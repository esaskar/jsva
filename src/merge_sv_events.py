#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
"""
Postprocess DELLY simple event calls to complex event calls by 
first merging simple calls across samples and then chaining
breakpoint clusters.

Orientation notation:
- Forward strand = + = 3' breakend    ----|
- Reverse strand = - = 5' breakend    |----

Author: Esa Pitkänen (esa.pitkanen@helsinki.fi)
"""

import sys, os, argparse, bisect, datetime, math, random
from collections import Counter
import numpy as np
import scipy
import circos

debug = 0

SV_DEL = 1
SV_DUP = 2
SV_INV = 3
SV_TR = 4

ORIENTATIONS = {"3to5" : 1,
                "5to3" : 2,
                "3to3" : 3,
                "5to5" : 4}

BREAKEND = {"3" : "+",
            "5" : "-"}

sample_ix = {}
ix_sample = {}
next_sample_ix = 1

def map_sample(sample):
    global sample_ix, ix_sample, next_sample_ix
    if sample not in sample_ix:
        sample_ix[sample] = next_sample_ix
        ix_sample[next_sample_ix] = sample
        next_sample_ix += 1
    return sample_ix[sample]

class Breakpoint:
    def __init__(self, chrom, pos, svtype, breakend):
        self.chrom = chrom
        self.pos = pos
        self.svtype = svtype
        self.ev = None
        self.breakend = breakend
    def get_mate(self):
        if self == self.ev.bp1:
            return self.ev.bp2
        else:
            return self.ev.bp1
    def same_pos_sample(self, o):
        "Return True iff two breakpoints have the same position in the same sample"
        return self.ev.sample == o.ev.sample and self.chrom == o.chrom and self.pos == o.pos
    def __cmp__(self, o):
        if o == None:
            return 1
        elif self.chrom == o.chrom:
            return cmp(self.pos, o.pos)
        else:
            return cmp(self.chrom, o.chrom)
    def __str__(self):
        return "bp(%s:%s:%s)" % (self.chrom, self.pos, self.breakend)
    def __repr__(self):
        return "bp(%s:%s:%s)" % (self.chrom, self.pos, self.breakend)
    def __hash__(self):
        return hash("%s:%s" % (self.chrom, self.pos))

class ComplexSVEvent:
    pass

class SimpleSVEvent:
    COLOR = {SV_DEL : "red",
             SV_DUP : "green",
             SV_INV : "blue",
             SV_TR : "orange"}
    EVENT_NAME = {SV_DEL : "del",
                  SV_DUP : "dup",
                  SV_INV : "inv",
                  SV_TR : "tr"}
    BP_PAIR_ORIENTATIONS = {SV_DEL : [("+", "-")],
                            SV_DUP : [("-", "+")],
                            SV_INV : [("-", "-"), ("+", "+")],
                            SV_TR : [("-", "-"), ("+", "+"), ("-", "+"), ("+", "-")]}
    def __init__(self, bp1, bp2, reads, qual, svid, svtype, sample):
        self.bp1 = bp1
        self.bp2 = bp2
        self.reads = reads
        self.qual = qual
        self.svid = svid
        self.svtype = svtype
        self.sample = sample
    def color(self):
        return SimpleSVEvent.COLOR[self.svtype]
    def __eq__(self, o):
        if o == None:
            return False
        return self.bp1 == o.bp1 and self.bp2 == o.bp2 and self.svtype == o.svtype
    def __str__(self):
        return "%s:%s:%s:%s(%s,%s)" % (SimpleSVEvent.EVENT_NAME[self.svtype], ix_sample[self.sample], self.reads, int(self.qual), self.bp1, self.bp2)
    def __repr__(self):
        return "%s(%s,%s)" % (SimpleSVEvent.EVENT_NAME[self.svtype], self.bp1, self.bp2)
    def __hash__(self):
        return hash(hash(self.bp1) + hash(self.bp2) + hash(self.svtype) + hash(self.sample))

class MergedSVEvent:
    "Collection of SimpleSVEvents merged due to SV type and breakpoint orientation identity, and breakpoint proximity."
    def __init__(self, events = []):
        self.events = set(events)
    def chrs_breakends(self):
        evs = list(self.events)
        return evs[0].bp1.chrom, evs[0].bp1.breakend, evs[0].bp2.chrom, evs[0].bp2.breakend
    def __str__(self):
        return "MergedSVEvent: %s" % (",".join(map(str, self.events)))
    def __hash__(self):
        evs = list(self.events)
        evs.sort()
        return hash(sum(map(hash, evs)))
    def __eq__(self, o):
        if o == None:
            return False
        else:
            n1 = len(self.events)
            n2 = len(o.events)
            if n1 != n2:
                return False
            else:
                return self.events == o.events
    def __cmp__(self, o):
        if o == None:
            return 1
        e1 = list(self.events)
        e2 = list(o.events)
        bp1 = min(e1[0].bp1, e1[0].bp2)
        bp2 = min(e2[0].bp1, e2[0].bp2)
        return cmp(bp1, bp2)

def map_breakend(svid, svtype):
    if svtype == SV_DEL:
        return "+", "-"
    elif svtype == SV_DUP:
        return "-", "+"
    elif svtype == SV_INV:
        if svid.split("_")[1] == "0":
            return "+", "+"
        else:
            return "-", "-"
    elif svtype == SV_TR:
        first, second = svid.split("_")[1].split("to")
        return BREAKEND[first], BREAKEND[second]
    else:
        print svid, svtype
        assert(0) # unknown svtype

def canonical_breakpoint_order(chr1, pos1, chr2, pos2):
    "Ensure that bp1 <= bp2. This should never happen with DELLY data."
    if chr1 > chr2 or (chr1 == chr2 and pos1 > pos2):
        sys.stderr.write("Warning: Swapping breakpoints: %s" ",".join([chr1, pos1, chr2, pos2]))
        return chr2, pos2, chr1, pos1
    else:
        return chr1, pos1, chr2, pos2

def chr_canonical(x):
    if x.startswith("chr"):
        return x[3:]
    else:
        return x

def create_event(vals, svtype, sample, args):
    if svtype == SV_TR:
        chr1, pos1, chr2, pos2, reads, qual, svid = vals
        chr1, chr2 = chr_canonical(chr1), chr_canonical(chr2)
        pos1, pos2 = int(pos1), int(pos2)
        chr1, pos1, chr2, pos2 = canonical_breakpoint_order(chr1, pos1, chr2, pos2)
        breakend_first, breakend_second = map_breakend(svid, svtype)
        bp1 = Breakpoint(chr1, pos1, svtype, breakend_first)
        bp2 = Breakpoint(chr2, pos2, svtype, breakend_second)
    else:
        chr1, pos1, end1, length, reads, qual, svid = vals
        chr1 = chr_canonical(chr1)
        pos1, end1 = int(pos1), int(end1)
        chr1, pos1, chr1, end1 = canonical_breakpoint_order(chr1, pos1, chr1, end1)
        breakend_first, breakend_second = map_breakend(svid, svtype)
        bp1 = Breakpoint(chr1, pos1, svtype, breakend_first)
        bp2 = Breakpoint(chr1, end1, svtype, breakend_second)
    if args.chrom != None and chr1 not in args.chrom:
        return None
    reads = int(reads)
    qual = float(qual)
    if reads < args.reads or qual < args.qual or qual > args.max_qual:
        return None
    else:
        ev = SimpleSVEvent(bp1, bp2, reads, qual, svid, svtype, sample)
        bp1.ev = ev
        bp2.ev = ev
        return ev

def add_breakpoint(bp, breakpoints):
    if bp.chrom not in breakpoints:
        breakpoints[bp.chrom] = []
    bisect.insort_left(breakpoints[bp.chrom], bp)

def add_breakpoints(event, breakpoints):
    add_breakpoint(event.bp1, breakpoints)
    add_breakpoint(event.bp2, breakpoints)

EDGE_EVENT = 1
EDGE_PROXIMAL = 2

def create_breakpoint_graph(events, breakpoints, args):
    G = {}
    n_event_edges = n_proximal_edges = 0
    for ev in events:
        G[ev.bp1] = {}
        G[ev.bp2] = {}
        G[ev.bp1][ev.bp2] = (EDGE_EVENT, 0)
        G[ev.bp2][ev.bp1] = (EDGE_EVENT, 0)
        n_event_edges += 1
    for chrom in breakpoints:
        prev_bpi = 0
        prev_bp = breakpoints[chrom][prev_bpi]
        for i, bp in enumerate(breakpoints[chrom]):
            while prev_bp != bp and bp.pos - prev_bp.pos > args.max_bp_dist:
                prev_bpi += 1
                prev_bp = breakpoints[chrom][prev_bpi]
            for j in range(prev_bpi, i):
                bpj = breakpoints[chrom][j]
                G[bpj][bp] = (EDGE_PROXIMAL, bp.pos - bpj.pos)
                G[bp][bpj] = (EDGE_PROXIMAL, bp.pos - bpj.pos)
                n_proximal_edges += 1

    print "%d event edges, %d proximal edges" % (n_event_edges, n_proximal_edges)

SV_FILE_SUFFIXES = {SV_DEL : "delly_pe.txt",
                    SV_INV : "invy_pe.txt",
                    SV_DUP : "duppy_pe.txt",
                    SV_TR : "jumpy_pe.txt"}

def get_sv_files(ddir, svtype):
    return ["%s/%s" % (ddir, fn) for fn in os.listdir(ddir) if os.path.isfile("%s/%s" % (ddir, fn)) and fn.endswith(SV_FILE_SUFFIXES[svtype])]

def parse_calls(ddir, svtype, events, breakpoints, sample, args):
    if SimpleSVEvent.EVENT_NAME[svtype] not in args.sv_type:
        return 0
    fns = get_sv_files(ddir, svtype)
#   print ddir, svtype, fns
    n_events = n_reads = n_pass = 0
    parsed = 0
    for fn in fns:
#        print "Opening", fn
        f = open(fn)
        if args.no_parse:
            continue
        for s in f:
            if s.startswith("#"):
                continue
            if args.max_events != None and n_events > args.max_events:
                break
            if ">" not in s:
                n_reads += 1
                continue
            n_events += 1
            ev = create_event(s.strip().split("\t"), svtype, sample, args)
            if ev != None:  # read and quality threshold check, chromosome check
                events.append(ev)
                add_breakpoints(ev, breakpoints)
                n_pass += 1
        parsed += 1
        f.close()
    if n_events == 0:
        fr = "N/A"
    else:
        fr = "%.2f%%" % (100.0 * n_pass / n_events)
    print "%s: %d/%d events (%s) accepted (%d total reads, %d files)" % (SimpleSVEvent.EVENT_NAME[svtype], n_pass, n_events, fr, n_reads, parsed)
    return n_events

def plot_circos(samples, events, breakpoints, args):
    for sample in samples:
        links = []
        histograms = []
        for ev in events:
            if ev.svtype == SV_TR or abs(ev.bp1.pos - ev.bp2.pos) > args.local_sv_len:
                links.append((ev.bp1.chrom, ev.bp1.pos, ev.bp2.chrom, ev.bp2.pos, ev.color()))
            else:
                histograms.append((ev.bp1.chrom, ev.bp1.pos, ev.bp2.pos, 1.0, ev.color()))
        ofn = "%s/%s_sv" % (args.output, sample)
        circos.plot(ofn, links, histograms)

class Sample:
    def __init__(self, sample, individual, tumortype, istumor, ddir):
        self.sample, self.individual, self.tumortype, self.istumor, self.ddir = sample, individual, tumortype, istumor, ddir

def parse_config(fn, max_samples):
    f = open(fn)
    samples = {}
    errors = 0
    c = 0
    for s in f:
        if s.startswith("#"):
            continue
        sample, individual, tumortype, istumor, ddir = s.strip().split("\t")
        if sample in samples:
            sys.stderr.write("%s: duplicate\n" % (sample))
            errors += 1
            continue
        if os.path.exists(ddir) == False:
            sys.stderr.write("%s: data dir not found\n" % (ddir))
            errors += 1
            continue
        if os.path.exists("%s/wgspipeline/delly" % (ddir)) == False and os.path.exists("%s/delly" % (ddir)) == False:
            sys.stderr.write("%s: delly dir not found\n" % (ddir))
            errors += 1
            continue
 
        samples[sample] = Sample(sample, individual, tumortype, istumor, ddir)
        c += 1
        if max_samples != None and c == max_samples:
            break

    if errors > 0:
        print "%d errors" % (errors)
        sys.exit(1)
    return samples

def find_nearby_breakpoints(bp, breakpoints, args):
    i = j = bisect.bisect_left(breakpoints, bp)
#    j = bisect.bisect_right(breakpoints, bp)
#    print "before:", i, j, bp, breakpoints[i], breakpoints[j]
    while i > 0 and bp.pos - breakpoints[i - 1].pos < args.merge_dist:
        i -= 1
    while j < len(breakpoints) - 2 and breakpoints[j + 1].pos - bp.pos <= args.merge_dist:
        j += 1
#    print "after:", i, j, bp, breakpoints[i], breakpoints[j]
    return breakpoints[i:j + 1]

def incompatible_breakpoints(ev1, ev2):
    global debug
    bp1a, bp1b = ev1.bp1, ev1.bp2
    bp2a, bp2b = ev2.bp1, ev2.bp2
    if debug: print "check incompatible-bp:", ev1, ev2
    # disregard ev1 == ev2
    if (bp1a.same_pos_sample(bp2a) and bp1b.same_pos_sample(bp2b)) or (bp1a.same_pos_sample(bp2b) and bp1b.same_pos_sample(bp2a)):
        if debug: print "same sample - returning"
        return True
    if debug: print ev1, bp1a, bp1b, ev2, bp2a, bp2b
    if bp1a.svtype != bp2a.svtype:
        if debug: print "different svtype - return"
        return True
    assert(bp1a.svtype == bp1b.svtype)
    assert(bp2a.svtype == bp2b.svtype)
    if not (bp1a.breakend == bp2a.breakend and bp1b.breakend == bp2b.breakend):
        if debug: print "different orientation - return"
        return True
    if debug: print "compatible bp"
    return False

def find_nearby_events(ev, events, breakpoints, args):
    global debug
    near_bp1 = find_nearby_breakpoints(ev.bp1, breakpoints[ev.bp1.chrom], args)
#    near_bp2 = find_nearby_breakpoints(ev.bp2, breakpoints[ev.bp2.chrom], args)
#    print ev.bp1, ev.bp2
    if debug: print "f_n_e:", ev.bp1, ev.bp2, "near_bp1:", near_bp1
    nearby_events = set()
#    print "ev.bp1", ev.bp1, "ev.bp2", ev.bp2
    for nbp1 in near_bp1:
        if debug: print "process", nbp1, nbp1.svtype, ev.svtype, ev.bp1.breakend
        if incompatible_breakpoints(ev, nbp1.ev):
            continue
        # nbp1 is near bp1, check if nbp2 is near bp2 
        if ev.bp2.chrom != nbp1.ev.bp2.chrom or abs(ev.bp2.pos - nbp1.ev.bp2.pos) > args.merge_dist:
            if debug: print "bps far - skip"
            continue
        # both breakpoints are nearby this event's breakpoints
        if debug: print "add nearby event:", nbp1.ev
        nearby_events.add(nbp1.ev)
    return nearby_events

def merged_event_breakpoint_stats(mev):
    bp1d, bp2d = [], []
    bend1 = bend2 = None
    reads = []
    quals = []
    for ev in mev.events:
        bp1d.append(ev.bp1.pos)
        bp2d.append(ev.bp2.pos)
        reads.append(ev.reads)
        quals.append(ev.qual)
        bend1 = ev.bp1.breakend
        bend2 = ev.bp2.breakend
    bp1d = np.array(bp1d)
    bp2d = np.array(bp2d)
    if bend1 == "+":
        bp1limit = scipy.amin(bp1d)
    else:
        bp1limit = scipy.amax(bp1d)
    if bend2 == "+":
        bp2limit = scipy.amin(bp2d)
    else:
        bp2limit = scipy.amax(bp2d)
    reads_median = int(scipy.median(reads))
    qual_median = int(scipy.median(quals))
    return int(bp1limit), int(bp2limit), int(bp2limit - bp1limit), scipy.mean(bp1d), scipy.amax(bp1d) - scipy.amin(bp1d), scipy.std(bp1d), scipy.mean(bp2d), scipy.amax(bp2d) - scipy.amin(bp2d), scipy.std(bp2d), reads_median, qual_median

def count_events(merged_events, args):
#    ev_count = Counter()
    o = open("%s/event_counts.txt" % (args.output), "w")
    o2 = open("%s/breakpoint_jitter.txt" % (args.output), "w")
    for me in merged_events:
#        print me
        bp1l, bp2l, svlen, m1, d1, sd1, m2, d2, sd2, reads_median, qual_median = merged_event_breakpoint_stats(me)
        o2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (me, m1, d1, sd1, m2, d2, sd2, reads_median, qual_median))
#    for ev in ev_count:
#        o.write("%s\t%d\n" % (ev, ev_count[ev]))

def merge_events(events, breakpoints, args):
    global debug
    merged_events = set()
    for ev in events:
        if debug: print ev
        # find events of same type whose both breakpoints are "close" this event's breakpoints
        near_events = find_nearby_events(ev, events, breakpoints, args)
        if debug: print "near_events: ", near_events
        near_events.add(ev) 
        merged_events.add(MergedSVEvent(near_events))  # store only unique merged events
    return merged_events

def find_breakpoint_clusters(bp_list, max_dist):
    i = 0
    n = len(bp_list)
    clusters = []
    cluster = []
    while i < n:
        if i > 0 and bp_list[i].pos - bp_list[i - 1].pos > max_dist:
            clusters.append(cluster)
            cluster = []
        cluster.append(bp_list[i])
        i += 1
    clusters.append(cluster)
    return clusters

def partition_breakpoints(bps, max_dist):
    chroms = {}
    # arrange mate breakpoints in chromosome-specific lists
    for bp in bps:
        bp2 = bp.get_mate()
        if bp2.chrom not in chroms:
            chroms[bp2.chrom] = [bp2]
        else:
            bisect.insort_left(chroms[bp2.chrom], bp2)
    # find clustered mate breakpoints
    partitions = []
    for chrom in chroms:
        partitions.extend(find_breakpoint_clusters(chroms[chrom], max_dist))
    return partitions

def merge_events_sweep_chrom(svtype, bp1or, bp2or, events, chrom_breakpoints, args):
    """
svtype -- merge events with this svtype
bp1or, bp2or -- require that BP1 and BP2 match these orientations
events -- list of simple SV events
chrom_breakpoints -- list of breakpoints to process
"""
    global debug
    merged_events = set()
    i = 0
    n = len(chrom_breakpoints)
#    print "bps:", chrom_breakpoints
    # i: first bp in event
    # j: current bp (sweep position)
    merged_breakpoints = set()
    while i < n:
#        print "start, i:", i
        bp1 = chrom_breakpoints[i]
        if bp1.svtype != svtype or bp1.breakend != bp1or or bp1.get_mate().breakend != bp2or or bp1 in merged_breakpoints:
            i += 1
            continue
        # i is now the first bp in the new event
        event_bps = [i]
#        breakpoint_mates = set([bp1.get_mate()])
        prev_bp = bp1
        j = i + 1
        # find the last breakpoint of the event
        while j < n:
            bp2 = chrom_breakpoints[j]
            # terminate if interval too long
            if bp2.pos - prev_bp.pos > args.merge_dist:         
                if debug: print "interval too long:", prev_bp.pos, "->", bp2.pos
                break
#            if bp2 in breakpoint_mates:
#                if debug: print "mate found:", bp2
#                break
            # skip breakpoints of different type or orientation
            if bp2.svtype != svtype or bp2.breakend != bp1or or bp1.get_mate().breakend != bp2or:   
                j += 1
                continue
            prev_bp = bp2
            event_bps.append(j)
#            breakpoint_mates.add(bp2.get_mate())
            j += 1
        # event_bps now contains clustered breakpoints
        # these breakpoints are next grouped into one or more events where both mate breakpoints are clustered
        if debug: print "event_bps:", map(lambda x: chrom_breakpoints[x].pos, event_bps)

        bp_partitions = partition_breakpoints(map(lambda x: chrom_breakpoints[x], event_bps), args.merge_dist)

        # create a MergedSVEvent for each mate breakpoint partition
        for p in bp_partitions:
            if debug: print "Partition:", p
            for bp in p:
                merged_breakpoints.add(bp)
                merged_breakpoints.add(bp.get_mate())
            mev = MergedSVEvent(map(lambda x: x.ev, p))
            if debug: print mev
            merged_events.add(mev)

        # continue sweep
        i = j
        if debug: print "i =", i

    return merged_events

def merge_events_sweep(events, breakpoints, args):
    global debug
    merged_events = set()
    for chrom in breakpoints:
        for svtype in SimpleSVEvent.EVENT_NAME:
            for bp1or, bp2or in SimpleSVEvent.BP_PAIR_ORIENTATIONS[svtype]:
                if debug: print "===== PROCESS", SimpleSVEvent.EVENT_NAME[svtype], bp1or, bp2or
                merged_events.update(merge_events_sweep_chrom(svtype, bp1or, bp2or, events, breakpoints[chrom], args))
    return merged_events

def get_shared_breakpoint_position(events, left = True):
    bps = []
    orient = None
    for x in events:
        if left:
            p = x.bp1
        else:
            p = x.bp2
        bps.append(p)
    if bps[0].breakend == "+":
        return reduce(lambda x, y: min(x, y), map(lambda x: x.pos, bps))
    else:
        return reduce(lambda x, y: max(x, y), map(lambda x: x.pos, bps))

class IndividualVariant:
    def __init__(self, individual):
        self.individual = individual
        self.normal_events = []  # list of normal events for this variant
        self.tumor_events = []   # list of tumor events for this variant
    def __cmp__(self, o):
        if o == None:
            return 1
        else:
            return cmp(self.individual, o.individual)
    def __str__(self):
        n = len(self.normal_events)
        t = len(self.tumor_events)
        tag = ""
        if n > 0:
            nquals = map(lambda x: x.qual, self.normal_events)
            nreads = map(lambda x: x.reads, self.normal_events)
            nqual = int(math.floor(1.0 * sum(nquals) / n))
            nread = int(math.floor(1.0 * sum(nreads) / n))
            nleftpos = get_shared_breakpoint_position(self.normal_events, left = True)
            nrightpos = get_shared_breakpoint_position(self.normal_events, left = False)
            nleftright = "%d-%d" % (nleftpos, nrightpos)
        else:
            nqual = nread = nleftright = "-"
        if t > 0:
            tquals = map(lambda x: x.qual, self.tumor_events)
            treads = map(lambda x: x.reads, self.tumor_events)
            tqual = int(math.floor(1.0 * sum(tquals) / t))
            tread = int(math.floor(1.0 * sum(treads) / t))
            tleftpos = get_shared_breakpoint_position(self.tumor_events, left = True)
            trightpos = get_shared_breakpoint_position(self.tumor_events, left = False)
            tleftright = "%d-%d" % (tleftpos, trightpos)
        else:
            tqual = tread = tleftright = "-"
        
        if n > 0:
            if t > 0:
                tag = "G" # germline
            else:
                tag = "N" # normal
        else:
            if t > 0:
                tag = "S" # somatic
            else:
                tag = "-" # No variant
        if tag != "-":
            if nread == "-":
                ntag = "-"
            else:
                ntag = "%s:%s:%s" % (nread, nqual, nleftright)
            if tread == "-":
                ttag = "-"
            else:
                ttag = "%s:%s:%s" % (tread, tqual, tleftright)
            tag = "%s/%s/%s" % (tag, ntag, ttag)
        return tag

def tumor_group_str(tumor_groups):
    #tumor_groups[sample.tumortype][sample.individual][sample.tumortype] += 1
    k = tumor_groups.keys()
    k.sort()
    v = []
    total_g = total_s = 0
    for tg in k:
        g = n = s = x = 0
        for indv in tumor_groups[tg]:
            nv = tumor_groups[tg][indv]["N"]
            tv = tumor_groups[tg][indv]["T"]
            if nv > 0:
                if tv > 0:
                    g += 1
                else:
                    n += 1
            else:
                if tv > 0:
                    s += 1
                else:
                    x += 1
        assert(x == 0)
#        print g, n, s, tumor_groups[tg]
        v.append("%d\t%d" % (g + n, s))
        total_g += g + n
        total_s += s

    return "%d\t%d\t%s" % (total_g, total_s, "\t".join(v))

def resolve_simple_events(merged_events, samples, args):
    inames = set()
    tgnames = set()
    ttindv = {}
    for six in samples:
        inames.add(samples[six].individual)
        tgnames.add(samples[six].tumortype)
        if samples[six].tumortype not in ttindv:
            ttindv[samples[six].tumortype] = set()
        ttindv[samples[six].tumortype].add(samples[six].individual)
    n_indv = len(inames)
    tgnames = list(tgnames)
    tgnames.sort()
    tgnames = "\t".join(map(lambda x: "%s_G/%d\t%s_S/%d" % (x, len(ttindv[x]), x, len(ttindv[x])), tgnames)) # create tumor group header: "CRC_G CRC_S myoma_G myoma_S ..."
    inames = list(inames)
    inames.sort()
    o = open("%s/simple_events.csv" % (args.output), "w")
    ho = open("%s/simple_events.header.txt" % (args.output), "w")
    hdr = "Chrom\tBP1\tBP2\tLen\tSVID\tSVType\tBP1Dir\tBP2Dir\tReadsMedian\tQualMedian\tGermline/%d\tSomatic/%d\t%s\t%s\tLengthMean\tBP1Chr\tBP1Mean\tBP1Diff\tBP1SD\tBP2Chr\tBP2Mean\tBP2Diff\tBP2SD\n" % (n_indv, n_indv, tgnames, "\t".join(inames))
    ho.write(hdr)
    ho.close()
    o.write("#Output of \"%s\" on %s\n" % (" ".join(sys.argv), datetime.datetime.now()))
    o.write(hdr)
    svid = 1
    for mev in merged_events:
        individual_variant = {}
        svtype = None
        tumor_groups = {}
        for six in samples:
            sample = samples[six]
            individual_variant[sample.individual] = IndividualVariant(sample.individual)
            if samples[six].tumortype not in tumor_groups:
                tumor_groups[samples[six].tumortype] = {}
        for ev in mev.events:
            assert(svtype == None or svtype == ev.svtype)
            svtype = ev.svtype
            sample = samples[ix_sample[ev.sample]]
            assert(sample.istumor == "N" or sample.istumor == "T")
            if sample.istumor == "N":
                individual_variant[sample.individual].normal_events.append(ev)
            else:
                individual_variant[sample.individual].tumor_events.append(ev)

            # count tumor types
            if sample.individual not in tumor_groups[sample.tumortype]:
                tumor_groups[sample.tumortype][sample.individual] = {"N" : 0, "T" : 0}
            tumor_groups[sample.tumortype][sample.individual][sample.istumor] += 1

        tgs = tumor_group_str(tumor_groups)

        ivs = individual_variant.values()
        ivs.sort()
        vstr = "\t".join(map(str, ivs))
        bp1l, bp2l, svlen, m1, d1, sd1, m2, d2, sd2, reads_median, qual_median = merged_event_breakpoint_stats(mev)
        chr1, or1, chr2, or2 = mev.chrs_breakends()
        if svtype != SV_TR:
            lenmean = m2 - m1
            assert(chr1 == chr2)
            chrom = chr1
        else:
            if chr1 != chr2:
                svlen = "-"
                lenmean = "-"
            else:
                lenmean = m2 - m1
            chrom = "%s/%s" % (chr1, chr2)
        svid_str = "SV%d" % (svid)
        o.write("%s\n" % ("\t".join(map(str, [chrom, bp1l, bp2l, svlen, svid_str,
                                              SimpleSVEvent.EVENT_NAME[svtype], 
                                              or1, or2,
                                              reads_median, qual_median,
                                              tgs, vstr, lenmean,
                                              chr1, m1, d1, sd1, 
                                              chr2, m2, d2, sd2]))))
        svid += 1

def main(args):
    try:
        os.makedirs(args.output)
    except:
        pass

    events = []
    breakpoints = {} # chr -> [bp]

    samples = parse_config(args.input, args.max_samples)
    zero_samples = []
    k = samples.keys()
    random.shuffle(k)
#    k = ["c919_1N", "c919_1T"]
    for i, sampleid in enumerate(k):
        sample = samples[sampleid]
        print "Parsing %d/%d:" % (i + 1, len(samples)), sample.sample, sample.ddir
        if os.path.exists("%s/wgspipeline/delly" % (sample.ddir)) == False:
            dellydir = "%s/delly" % (sample.ddir)
        else:
            dellydir = "%s/wgspipeline/delly" % (sample.ddir)
        sid = map_sample(sample.sample)
        n = 0
        n += parse_calls(dellydir, SV_DEL, events, breakpoints, sid, args)
        n += parse_calls(dellydir, SV_INV, events, breakpoints, sid, args)
        n += parse_calls(dellydir, SV_DUP, events, breakpoints, sid, args)
        n += parse_calls(dellydir, SV_TR, events, breakpoints, sid, args)
        if n == 0:
            zero_samples.append(sampleid)
        print "%d raw events, %d filtered events in all samples so far" % (n, len(events))
        sys.stdout.flush()
    ozf = open("zero_event_samples.txt", "w")
    ozf.write("\n".join(zero_samples))
    ozf.close()
    print "Merging SV events..."
#    merged_events = merge_events(events, breakpoints, args)
    merged_events = merge_events_sweep(events, breakpoints, args)
    merged_events = list(merged_events)
    merged_events.sort()
    print "Counting events..."
    count_events(merged_events, args)
    resolve_simple_events(merged_events, samples, args)
#    exit()
#    create_breakpoint_graph(events, breakpoints, args)
#    plot_circos(samples, events, breakpoints, args)
#    print breakpoints

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", help = "Input config file", required = True)
    p.add_argument("-o", "--output", help = "Output dir", required = True)
    p.add_argument("-r", "--reads", help = "Minimum read support [%(default)s]", default = 3, type = int)
    p.add_argument("-q", "--qual", help = "Minimum quality [%(default)s]", default = 10, type = float)
    p.add_argument("--max-qual", help = "Maximum quality (>253 are bad) [%(default)s]", default = 253, type = float)
    p.add_argument("-d", "--max-bp-dist", help = "Maximum breakpoint distance for joining two breakpoints in breakpoint graph [%(default)s]", default = 200, type = int)
    p.add_argument("-l", "--local-len", help = "Minimum SV length to show as links in circos [%(default)s]", default = 10000, dest = "local_sv_len", type = int)
    p.add_argument("-m", "--merge-dist", help = "Maximum distance to merge SV events [%(default)s]", default = 500, type = int)
    p.add_argument("-x", "--max-events", help = "Debug: max number of events", default = None, type = int)
    p.add_argument("--max-samples", help = "Maximum samples to parse [unlimited]", default = None, type = int)
    p.add_argument("--sv-type", help = "SV types to process, comma-separated [all]", default = None, type = str)
    p.add_argument("--chrom", help = "List of chromosomes to process, comma-separated [all]", default = None, type = str)
    p.add_argument("--debug", action = "store_true")
    p.add_argument("--no-parse", action = "store_true")
    args = p.parse_args()
    if args.debug:
        debug = 1
    if args.sv_type != None:
        args.sv_type = args.sv_type.split(",")
        sys.stderr.write("SV types to parse: %s\n" % (",".join(args.sv_type)))
    else:
        args.sv_type = SimpleSVEvent.EVENT_NAME.values()
    if args.chrom != None:
        args.chrom = args.chrom.split(",")
        sys.stderr.write("Parsing chromosomes %s\n" % (",".join(args.chrom)))
    if args.max_events != None:
        sys.stderr.write("Warning: parsing maximum of %d events / sample / svtype\n" % (args.max_events))
    main(args)
