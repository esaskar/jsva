#!/usr/bin/env python

import sys, argparse, os, gzip, errno
import circos
from unravel import SampleBreakpoint

SVTYPE_COLOR = {"del" : "red", "inv" : "blue", "dup" : "green", "tr" : "orange"}

def conv_reads(x):
    x = int(x)
    if x > 10:
        x = 10
    return x

def closer_bp(bp, bp1, bp2):
    if bp[0] == bp1[0]:
        if bp[0] != bp2[0]:
            return bp1
        elif abs(bp[1] - bp1[1]) < abs(bp[1] - bp2[1]):
            return bp1
        else:
            return bp2
    elif bp[0] == bp2[0]:
        return bp2
    else:
        return None

LINK_COLORS = {20 : "red", 
               12 : "dorange",
               9 : "orange",
               7 : "green",
               5  : "blue",
               0  : "vlgrey"}

def reads_to_color(reads):
    k = LINK_COLORS.keys()
    k.sort(lambda x, y: cmp(y, x))
    for i in k:
        if reads >= i:
            return LINK_COLORS[i]
    return "black"

def reads_to_text_size(reads):
    if reads >= 9:
        return "26p"
    elif reads >= 5:
        return "22p"
    else:
        return "18p"

def main(args):
    if args.input == None or args.input == "-":
        f = sys.stdin
    else:
        f = open(args.input)
    sample_bps = {}
    while 1:
        s = f.readline()
        if s.startswith("#"):
            continue
        else:
            hdr = s.strip()
            break

    for s in f:
        sbp = SampleBreakpoint(s)
        if sbp.sample not in sample_bps:
            sample_bps[sbp.sample] = []
        sample_bps[sbp.sample].append(sbp)

    text_labels = set()

    keys = sample_bps.keys()
    keys.sort()
#    keys = ["c1049_1"]

    # Plot individual samples
    for sample in keys:
        ofn = "%s/%s" % (args.output, sample)
        links = []
        histogram = []
        text = []
        for bp in sample_bps[sample]:
            # links: chr1, pos1, chr2, pos2, color
            # histogram: chr1, pos1, pos2, value, color
            if bp.svtype == "tr" or int(bp.length) > args.min_link_length:
                # "long" SV connection
                thickness = conv_reads(bp.treads)
                gene, genedist = bp.gene
                label = bp.get_label(args.gene_dist)
                reads = int(bp.treads)
                if reads >= args.min_reads_label and label != None and label not in text_labels:
                    text.append(circos.Annotation(bp.chrom, int(bp.start), int(bp.end), label, SVTYPE_COLOR[bp.svtype], reads_to_text_size(reads)))
                    text_labels.add(label)
                color = reads_to_color(int(bp.treads))
                links.append(circos.Link(bp.mbp1chrom, float(bp.mbp1pos), bp.mbp2chrom, float(bp.mbp2pos), color, thickness, label))
            else:
                # "short" SV connection
                reads = int(bp.treads)
                color = reads_to_color(int(bp.treads))
                norm_reads = 1.0 * min(reads, 6) / 6
                z = min(50, 50 - reads)
#                print bp.treads, norm_reads, color, z
                histogram.append(circos.Region(bp.chrom, float(bp.mbp1pos), float(bp.mbp2pos), value = norm_reads, color = color, z = z))
                gene, genedist = bp.gene
                label = bp.get_label(args.gene_dist)
                if reads >= args.min_reads_label and label != None and label not in text_labels:
                    pos = int(0.5 * (float(bp.mbp1pos) + float(bp.mbp2pos)))
                    text.append(circos.Annotation(bp.chrom, pos, pos, label = label, color = color, size = reads_to_text_size(int(bp.treads))))
                    text_labels.add(label)
        print "Plotting", ofn
        circos.plot(ofn, links, histogram, text)

#    print sample_bps

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input")
    p.add_argument("-o", "--output", required = True)
    p.add_argument("--min-link-length", default = 10000, type = int, help = "[%(default)s]")
    p.add_argument("--gene-dist", default = 10000, type = int, help = "[%(default)s]")
    p.add_argument("--min-reads-label", default = 5, type = int, help = "Minimum number of reads to annotate SV with gene name [%(default)s]")
    args = p.parse_args()
    main(args)
