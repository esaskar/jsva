#!/usr/bin/env python2.7
"""
chain.py - Chain adjacent breakpoints with compatible orientations to find complex structural events from jointly analysed DELLY data.

Event calls:
* DEL: Simple deletions (two breakpoints, A --> B <-- C)
* DUP: Tandem duplication / interchromosomal rearrangement with only two observed breakpoints A <-- B B --> C
* INV: Inversion (>=4 breakpoints)  A --> <== inv(B) --> <== C
* CTX: Interchromosomal translocation / complex event spanning at least two chromosomes
* ITX: Intrachromosomal rearrangement / inversion with only 2 observed breakpoints

Annotations:
* Complex: >4 breakpoints
* MultiSample: called in >1 sample (possible germline polymorphism/reference error)
* LowConf: median tumor reads <= --low-conf-reads
* Repeat: median breakpoint distance to RepeatMasker region < --dist-interpret-repeat

- 1k merge window @ merit:  462 s, 11.6 GB
- 10k merge window @ merit: 
- Takes ~11 GB memory on the 20140513b dataset
"""

import numpy as np # median
import sys, breakpoint, argparse, gzip, datetime

from breakpoint import Breakpoint, BreakpointSet

TAG_REPEAT = "Repeat"
TAG_COMPLEX = "Complex"
TAG_LOWCONF = "LowConf"
TAG_MULTISAMPLE = "MultiSample"

EVENT_CTX = "CTX"
EVENT_ITX = "ITX"
EVENT_DUP = Breakpoint.SV_TYPE_NAME[Breakpoint.SV_DUP].upper()
EVENT_DEL = Breakpoint.SV_TYPE_NAME[Breakpoint.SV_DEL].upper()
EVENT_INV = Breakpoint.SV_TYPE_NAME[Breakpoint.SV_INV].upper()

CHAINED_IN = 1      # breakpoints --| |--
CHAINED_OUT = 2     # breakpoints |-- --|
CHAINED_PAIR = 3    # breakpoints chained due to pairing 

CHAIN_TYPES = {CHAINED_IN : "IN",
               CHAINED_OUT : "OUT",
               CHAINED_PAIR : "PAIR"}

def openf(fn):
    try:
        return gzip.open(fn)
    except:
        return open(fn)

logf = sys.stderr

def logmsg(s):
    logf.write("%s: %s\n" % (datetime.datetime.now(), s))
    logf.flush()

def find_components(G):
    Q = map(lambda x: (x, None), G.keys())
    comp = {}
    next_comp = 1
    while len(Q) > 0:
        u, c = Q.pop()
        if u not in comp:
            if c != None:
                comp[u] = c
            else:
                comp[u] = next_comp
                next_comp +=1
            for v in G[u]:
                Q.append((v, comp[u]))
    return comp

def is_cyclic(G, comps):
    num_edges = {}
    num_nodes = {}
    for u in comps:
        c = comps[u]
        if c not in num_edges:
            num_edges[c] = len(G[u])
            num_nodes[c] = 1
        else:
            num_edges[c] += len(G[u])
            num_nodes[c] += 1
    cyclic = {}
    for c in num_nodes:
        nn = num_nodes[c]
        ne = num_edges[c] / 2
        cyclic[c] = ne >= nn
    return cyclic

def component_sizes(comp):
    comp_size = {}
    for u in comp:
        c = comp[u]
        if c not in comp_size:
            comp_size[c] = 1
        else:
            comp_size[c] += 1
    return comp_size

def chain(bps, outf, max_chain_dist, **kwargs):
#    d = kwargs["dist"]
    bp_buffer = []
    G = {}

    pairs = {}

    processed_chrs = set()
    
    for i, bp in enumerate(bps):
#        if i > 200000:
#            break
        if bp.chrom not in processed_chrs:
            processed_chrs.add(bp.chrom)
            logmsg("Chromosome %s..." % (bp.chrom))
        if bp.svid not in pairs:
            pairs[bp.svid] = [bp]
        else:
            pairs[bp.svid].append(bp)
        while len(bp_buffer) > 0 and (bp.pos - bp_buffer[0].pos > max_chain_dist or bp.chrom != bp_buffer[0].chrom):
            bp_buffer.pop(0)
#        print "Processing", bp
        for bp2 in bp_buffer:
            # check whether bp2-bp should be chained
#            print "Check", bp2
            # bp2.pos < bp.pos
            n_common = len(bp.somatic_samples.intersection(bp2.somatic_samples))
            if n_common == 0:
                continue
            if (bp2.orient != bp.orient):  # |--bp2-- --bp--|  or  --bp2--| |--bp--
#                print "Chain", bp2, bp, ",".join(bp.somatic_samples.intersection(bp2.somatic_samples))
                if bp not in G:
                    G[bp] = {}
                if bp2 not in G:
                    G[bp2] = {}

                if bp2.orient == Breakpoint.BREAKEND_3P:
                    chain_type = CHAINED_IN             # ---| |---
                else:
                    chain_type = CHAINED_OUT            # |--- ---|

                d = abs(bp.pos - bp2.pos)
                G[bp][bp2] = (chain_type, d)
                G[bp2][bp] = (chain_type, d)
        bp_buffer.append(bp)
#        print "buffer", len(bp_buffer)

    logmsg("Chaining paired breakpoints...")

    for svid in pairs:
        if len(pairs[svid]) != 2:
#            logmsg("Warning: breakpoints for %s: %d != 2" % (svid, len(pairs[svid])))
            pass
        else:
            bp1, bp2 = pairs[svid]
            if bp1 not in G:
                G[bp1] = {}
            if bp2 not in G:
                G[bp2] = {}
            if bp1.chrom == bp2.chrom:
                d = abs(bp1.pos - bp2.pos)
            else:
                d = "-"
            G[bp1][bp2] = (CHAINED_PAIR, d)
            G[bp2][bp1] = (CHAINED_PAIR, d)

    return G

def component_analysis(bps, G, outf, **kwargs):
    logmsg("Finding connected components...")
    comps = find_components(G)
    logmsg("Calculating component sizes...")
    cs = component_sizes(comps)
#    print cs
    logmsg("Finding cyclic components...")
    cyclic = is_cyclic(G, comps)
#    print cyclic

    component_bps = {}
    for bp in comps:
        comp = comps[bp]
        if comp not in component_bps:
            component_bps[comp] = []
        component_bps[comp].append(bp)

    outf.write("Event\tEventType\tFilter\tnBreakpoints\tIsCyclic\tnSamples\tUniqueSamples\tSVIDs\tChromosomes\tSVTypes\tMinTumorReads\tMedianTumorReads\tRepeats\tRegions\tEdges\n")
    keys = component_bps.keys()
    keys.sort()
    c = 0
    # each component is a structural event, output consists of one event/line 
    for comp in keys:
        if cs[comp] < kwargs["min_event_size"]:
            continue

        c += 1
        cbps = component_bps[comp]
        if cyclic[comp]:
            cyclicstr = "Yes"
        else:
            cyclicstr = "No"
        somatic_samples = set()
        chrs = set()
        svtypes = []
        treads = []
        svids = set()
        repeats = []
        chr_regions = {}  # chr -> (min, max)
#        print "EV%d: %s" % (c, cbps)
        comp_edges = []
        for bp in cbps:
            for bp2 in G[bp]:
                if bp2 in component_bps[comp]:
                    comp_edges.append("%s:%s-%s:%s/%s/%s" % (bp.chrom, bp.pos, bp2.chrom, bp2.pos, CHAIN_TYPES[G[bp][bp2][0]], G[bp][bp2][1]))

            somatic_samples.update(bp.somatic_samples)
            chrs.add(bp.chrom)
            svtypes.append(bp.svtype)
            treads.extend(bp.tumor_read_counts())
            svids.add(bp.svid)
            repeats.append((bp.data[bps.header_ix["RepeatMasker_FAMILY"]], bp.data[bps.header_ix["RepeatMasker_DIST"]]))
            if bp.chrom not in chr_regions:
                chr_regions[bp.chrom] = (bp.pos, bp.pos)
            else:
                region = chr_regions[bp.chrom]
                chr_regions[bp.chrom] = (min(region[0], bp.pos), max(region[1], bp.pos))
        somatic_samples = list(somatic_samples)
        somatic_samples.sort()
        sstr = ",".join(somatic_samples)
        chrs = list(chrs)
        chrs.sort()
        cstr = ",".join(chrs)
        svtstr = ",".join(map(Breakpoint.SV_TYPE_NAME.get, svtypes))
        treads.sort()
        svids = list(svids)
        svidstr = ",".join(svids)
        if len(treads) > 0:
            median_treads = np.median(treads)
        min_treads = min(treads)

        repeat_dists = map(lambda x: abs(int(x[1])), repeats)
        repeat_dists.sort()
        median_repeat_dist = np.median(repeat_dists)

        interpretation = "Unknown"
        svtypes = set(svtypes)
        if len(chrs) == 1:
            if len(svtypes) == 1:
                svtype = list(svtypes)[0]
                if svtype == Breakpoint.SV_TR:
                    interpretation = [EVENT_CTX]    # happens if breakpoints have been filtered prior to chaining, still mark as an interchromosomal event
                elif (svtype == Breakpoint.SV_INV and len(cbps) == 4) or (svtype == Breakpoint.SV_DUP and len(cbps) == 2) or svtype == Breakpoint.SV_DEL:
                    # require that inversions have four breakpoints and duplications two
                    interpretation = [Breakpoint.SV_TYPE_NAME[svtype].upper()]  # inv or dup
                else:
                    # everything else that happens within a single chromosome is an intrachromosomal event
                    interpretation = [EVENT_ITX]  # <-- <--, --> --> and <-- --> where more than 2 breakpoints are seen
            else:
                # events with multiple types of calls happening in one chromosome are just intrachromosomal events
                interpretation = [EVENT_ITX]
        else:
            # events involving multiple chromosomes are interchromosomal events
            interpretation = [EVENT_CTX]

        filter_str = []
        if median_repeat_dist < kwargs["dist_interpret_repeat"]:
            filter_str.append(TAG_REPEAT)

        if median_treads <= kwargs["low_conf_reads"]:
            filter_str.append(TAG_LOWCONF)

        if len(cbps) > 4:
            filter_str.append(TAG_COMPLEX)

        if len(somatic_samples) > 1:
            filter_str.append(TAG_MULTISAMPLE)

        interpretation = ",".join(interpretation)
        filter_str = ",".join(filter_str)

        region_str = []
        k = chr_regions.keys()
        k.sort()
        for chrom in k:
            region_str.append("%s:%d-%d" % (chrom, chr_regions[chrom][0], chr_regions[chrom][1]))
        region_str = ",".join(region_str)

        outf.write("EV%d\t%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n" % (c, interpretation, filter_str, len(cbps), cyclicstr, len(somatic_samples), sstr, svidstr, cstr, svtstr, min_treads, median_treads, ",".join(map(lambda x: "%s:%s" % (x[0], x[1]), repeats)), region_str, ",".join(comp_edges)))

def breakpoint_interdistances(bps, outf, **kwargs):
    last_pos = 0
    last_chrom = None
    for bp in bps:
        if bp.chrom != last_chrom:
            last_chrom = bp.chrom
            last_pos = bp.pos
        else:
            outf.write("%s\t%d\t%d\n" % (bp.chrom, bp.pos, bp.pos - last_pos))
            last_pos = bp.pos

def main(infn = None, outfn = None, **kwargs):
    inf = sys.stdin if infn == "-" else openf(infn)
    outf = sys.stdout if outfn == "-" else open(outfn, "w")
    bps = breakpoint.BreakpointSet(inf, max_germline_af = 0.0)
    if kwargs["interdistances"]:
        breakpoint_interdistances(bps, outf, **kwargs)
    else:
        logmsg("Chaining breakpoints...")
        G = chain(bps, outf, **kwargs)
        logmsg("Performing component analysis of the breakpoint chain...")
        component_analysis(bps, G, outf, **kwargs)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", dest = "infn", default = "-")
    p.add_argument("-o", "--output", dest = "outfn", default = "-")
    p.add_argument("-d", "--distance", help = "Maximum chaining distance [%(default)d]", dest = "max_chain_dist", metavar = "DIST", default = 1000, type = int)
    p.add_argument("--interdistances", action = "store_true")
    p.add_argument("--min-event-size", help = "Minimum event size in breakpoints [%(default)d]", default = 1, type = int)
    p.add_argument("--dist-interpret-repeat", help = "Maximum median distance between event breakpoints and closest RepeatMasker region to interpret event as repeat [%(default)d]", default = 100, type = int)
    p.add_argument("--low-conf-reads", help = "Maximum median tumor read count to call event low confidence [%(default)d]", default = 4, type = int)
    args = p.parse_args()
    main(**vars(args))
