#!/usr/bin/env python

import sys, argparse, os, gzip, errno
import circos_links_only as circos
from breakpoint import Breakpoint, BreakpointSet
import chain

def canonical_sample(s):
    if s.startswith("Fam_"):
        s = s[4:]
    v = s.strip().split("_")
    if len(v) == 1 or len(v[1]) > 3:
        # extraction number
        return "%s_1" % (v[0])
    else:
        return "%s_%s" % (v[0], v[1])

def parse_events(f, args):
    hdr = None
    sample_events = {}

    for s in f:
        if s.startswith("#"):
            continue
        if hdr == None:
            hdr = s.strip().split("\t")
            continue
        evid, evtype, filterstr, n_bp, iscyclic, n_samples, samples, svids, chroms, svtypes, min_tumor_reads, median_tumor_reads, repeats, regions, edges = s.strip().split("\t")
#        bps, chain_type, dist = edges.split("/")
#        bp1, bp2 = map(lambda x: x.split(":"), bps.split("-"))

        tags = filterstr.split(",")
        if args.allow_low_conf == False and chain.TAG_LOWCONF in tags:
            continue
        if args.allow_multisample == False and chain.TAG_MULTISAMPLE in tags:
            continue
        if args.complex_only and chain.TAG_COMPLEX not in tags:
            continue
        if args.disallow_repeats and chain.TAG_REPEAT in tags:
            continue

        samples = samples.split(",")
        for sample in samples:
            sample = canonical_sample(sample)
            if sample not in sample_events:
                sample_events[sample] = []
            sample_events[sample].append((evid, evtype, tags, min_tumor_reads, median_tumor_reads, edges))

    return sample_events


COLOR_COMPLEX = "152,78,163,0.3"

EVENT_COLORS = {chain.EVENT_ITX : "255,127,0,0.3",
                chain.EVENT_CTX : "255,127,0,0.3",
                chain.EVENT_DEL : "77,175,74,0.3",
                chain.EVENT_DUP : "228,26,28,0.3",
                chain.EVENT_INV : "55,126,184,0.3"}

#228,26,28   # red
#55,126,184  # blue
#77,175,74   # green
#152,78,163  # purple
#255,127,0   # orange

def treads_to_lw(treads):
    return 1 + min(max(1, int(0.7 * treads)), 20)

def draw_events(events, cnv, args):
    for sample in events:
#        sample = "c300_1"
        print sample
        links = []
        complex_chroms = set()
        for evid, evtype, tags, min_tumor_reads, median_tumor_reads, edges in events[sample]:
#            print sample, evid, edges
#            print edges
            edges = edges.split(",")
            iscomplex = chain.TAG_COMPLEX in tags
            if iscomplex:
                color = COLOR_COMPLEX
                z = 2
            else:
                color = EVENT_COLORS[evtype]
                z = 1
            thickness = treads_to_lw(int(median_tumor_reads))
            for edge in edges:
                bps, chain_type, dist = edge.split("/")
#                if iscomplex == False and chain_type != chain.CHAIN_TYPES[chain.CHAINED_PAIR]:
#                    continue
                if chain_type != chain.CHAIN_TYPES[chain.CHAINED_PAIR]:
                    continue
                bp1, bp2 = map(lambda x: x.split(":"), bps.split("-"))
#                print bp1, bp2, chain_type, dist
                chr1 = circos.conv_chr(bp1[0])
                chr2 = circos.conv_chr(bp2[0])

                if iscomplex:
                    complex_chroms.add(chr1)
                    complex_chroms.add(chr2)

                label = ""
                links.append(circos.Link(chr1, int(bp1[1]), chr2, int(bp2[1]), color, thickness, label, z))
 
        histogram = []
        if sample in cnv:
#            print "%d CNV regions" % (len(cnv[sample]))
            for chrom, start, end, aitype, baf, logrr in cnv[sample]:
                z = 1
                if aitype == "Loss":
                    color = "green"
                    value = logrr + 0.5
                elif aitype == "Gain":
                    color = "red"
                    value = logrr + 0.5
                elif aitype == "LOH":
                    r = g = (1.0 - baf) * 2
#                    print baf, r
                    value = logrr + 0.5
                    color = "(%d,%d,250)" % (int(255 * r), int(255 * g))
                    z = 3
                else:
                    color = "black"
                    value = 0.5
#                print sample, chrom, start, end, aitype, len(cnv[sample])
                histogram.append(circos.Region(circos.conv_chr(chrom), start, end, value, color, z = z))

        ofn = "%s/%s.png" % (args.odir, sample)
        circos.plot(ofn, links, histogram, [], chroms = None)  # plot all chroms

        if len(complex_chroms) > 0:
            ofn = "%s/%s.complex.png" % (args.odir, sample)
            circos.plot(ofn, links, histogram, [], chroms = complex_chroms)

#        sys.exit()

def cmp_region(x, y):
    if x[0] == y[0]:
        return cmp(x[1], y[1])
    else:
        return cmp(x[0], y[0])

def parse_cnv_snpchip(cnvfn):
    f = gzip.open(cnvfn)
    cnv = {}
    for s in f:
        chrom, start, end, aitype, baf, logrr, sample = s.strip().split("\t")
        ss = sample
        sample = canonical_sample(sample)
        if sample not in cnv:
            cnv[sample] = []
        cnv[sample].append((chrom, int(start), int(end), aitype, float(baf), float(logrr)))

#    for sample in cnv:
#        cnv[sample].sort(cmp_region)

    return cnv

def parse_cnv_varscan(cnvdir):
    sys.stderr.write("Parsing CNV data in VarScan format...\n")
    fns = os.listdir(cnvdir)
    cnv = {}
    for cnvfn in fns:
        if cnvfn.endswith(".bed.gz") == False:
            continue
        f = gzip.open("%s/%s" % (cnvdir, cnvfn))
        sample = canonical_sample(cnvfn)
        print cnvfn, sample
        if sample not in cnv:
            cnv[sample] = []
        for s in f:
            if s[0] == "#":
                continue
            chrom, start, end, numsegs, logcnratio, call = s.strip().split("\t")
            if call == "loss":
                call = "Loss"
            elif call == "gain":
                call = "Gain"
            elif call == "neutral":
                call = "LOH" # not true, we don't know AI status
            baf = 0.5 # incorrect
            cnv[sample].append((chrom, int(start), int(end), call, baf, float(logcnratio)))

#    for sample in cnv:
#        cnv[sample].sort(cmp_region)

    return cnv

    
def main(args):
    if args.input:
        if args.input.endswith(".gz"):
            f = gzip.open(args.input)
        else:
            f = open(args.input)
    else:
        f = sys.stdin
    if args.cnv:
        if args.varscan:
            cnv = parse_cnv_varscan(args.cnv)
        else:
            cnv = parse_cnv_snpchip(args.cnv)
    else:
        cnv = None
    events = parse_events(f, args)
    draw_events(events, cnv, args)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", help = "Input: chained SVs")
    p.add_argument("-o", "--odir", help = "Output dir for circos plots", required = True)
    p.add_argument("--allow-low-conf", help = "Allow low confidence events", action = "store_true")
    p.add_argument("--allow-multisample", help = "Allow multisample calls", action = "store_true")
    p.add_argument("--complex-only", help = "Show only complex calls", action = "store_true")
    p.add_argument("--disallow-repeats", help = "Do not show calls hitting repeats", action = "store_true")
    p.add_argument("--cnv", help = "Use CNV data from file")
    p.add_argument("--varscan", help = "CNV data in varscan format [array]", action = "store_true")
    args = p.parse_args()
    main(args)
