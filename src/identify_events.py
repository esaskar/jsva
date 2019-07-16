#!/usr/bin/env python2
"""
Identify simultaneous rearrangement events spanning multiple breakpoints in structural variation data.

This script builds an event graph from paired-end breakpoint and copy number variation data.
Two breakend nodes are joined in the graph iff
   1) breakends are paired mates in SV data or
   2) breakends are in close proximity and (optionally) have opposite orientations (head-tail or tail-head)

We assume that proximal breakends in (+,-) (--| |--) orientation have emerged from the same rearrangement event.
Likewise, proximal breakends in (-,+) (|-- --|) orientation are assumed to belong to the same fragment
and thus be involved in the same event.

Rearrangements events consist of breakends that have plausibly occurred at the same time.
Such events are identified by finding connected components of the event graph.
An event is called complex if it involves at least three breakpoints
(six breakends, or four breakends such that the rearrangement is not a balanced translocation).

Before constructing the event graph, deleted regions are identified in copy number data
and the breakend positions are adjusted to take into account deletions (distance
between any two breakend either decreases or remains the same).

A permutation test is performed to calculate a p-value for the result.
Positions of paired-end breakends are randomly shuffled within the original chromosome
while preserving mate pairing. Copy number data is not permuted.
Test statistic is by default the number of breakends associated with called complex events.

Author: Esa Pitkanen (esa.pitkanen@cs.helsinki.fi)
"""

import sys, argparse, tempfile, subprocess, datetime, re, StringIO, random, bisect, os, gzip, locale
import scipy.stats as stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import unravel  # to parse jsva result files

# Complex event test statistic
STATISTIC_EVENT_COUNT = 1               # Number of called complex events
STATISTIC_EVENT_BREAKPOINT_COUNT = 2    # Number of breakpoints in called complex events
DEFAULT_STATISTIC = STATISTIC_EVENT_BREAKPOINT_COUNT

# Default threshold for calling deletion in CNV data
DEFAULT_CNV_DEL_THRESHOLD = -0.4

# What to do when a breakpoint is inside deletion region
STRATEGY_DEL_BP_REPOSITION = 1   # reposition breakpoint to deletion region boundary
STRATEGY_DEL_BP_DELETE = 2       # delete breakpoint
DEFAULT_DEL_BP_STRATEGY = STRATEGY_DEL_BP_DELETE

# Distance threshold to join two breakpoints in the event graph
DEFAULT_JOIN_DISTANCE = 10000

# Minimum event cluster size to call complex event for statistics
DEFAULT_MIN_COMPLEX_SIZE = 5

# Default randomization rounds
DEFAULT_ROUNDS = 100

# Accepted input breakpoint formats
INPUT_ANNOTATED = 1          # "event chr pos T/H [comment]" or "deletion" or empty line
INPUT_BREAKDANCER_FR = 2     # "chr1 pos1 orientation1 chr2 pos2 orientation2", orientation=F/R
INPUT_BREAKDANCER_RAW = 3    # "chr1 pos1 orientation1 chr2 pos2 orientation2", orientation=X+Y-
INPUT_CG_RAW = 4             # Complete Genomics data
INPUT_DELLY = 5              # DELLY result file (JUMPY paired-end)
INPUT_DELLY_FILTERED = 6     # Concatenated and filtered DELLY result files
INPUT_JSVA = 7               # Joint structural variant analysis result (from unravel.py)
DEFAULT_INPUT_TYPE = INPUT_BREAKDANCER_FR

DEFAULT_EXCLUDE_REPEATS = ["centr"]
DEFAULT_CHROMS = map(str, range(1, 23)) + ["X", "Y", "MT"]

re_br_or = re.compile("(\d+)\+(\d+)\-")

LINK_END = 1
LINK_EVENT = 2

CNV_DATA_FORMAT_RAW = 1
CNV_DATA_FORMAT_CALLED = 2
CNV_DATA_FORMAT_VARSCAN = 3

BREAKEND_3P = "F"  # ---|
BREAKEND_5P = "R"  #      |---

JSVA_BREAKEND_ORIENTATION = {"+" : BREAKEND_3P, "-" : BREAKEND_5P}

class Node:
    def __init__(self, node_id, event_id, chrom, pos, end):
        self.node_id = node_id
        self.event_id = event_id
        self.chrom = chrom
        self.pos = pos
        self.end = end
        self.adj_pos = None # adjusted position after considering deleted regions in contract()
        self.reads = (None, None)   # number of supporting normal and tumor reads
        self.gene = (None, None)    # nearest gene and distance
        self.repeat = (None, None)  # nearest repeat and distance

    def __str__(self):
        return "%s %s %s %s->%s %s" % (self.node_id, self.event_id, self.chrom, self.pos, self.adj_pos, self.end)

    def __repr__(self):
        return "(%s %s %s %s %s)" % (self.node_id, self.event_id, self.chrom, self.pos, self.end)

class ContractNode:
    NODE_BP = 1
    NODE_CNV_START = 2
    NODE_CNV_END = 3
    def __repr__(self):
        if self.bptype == ContractNode.NODE_BP:
            s = "bp"
        elif self.bptype == ContractNode.NODE_CNV_START:
            s = "del_start"
        elif self.bptype == ContractNode.NODE_CNV_END:
            s = "del_end"
        else:
            assert(0) # invalid bptype
        return "(%s, %s)" % (s, self.pos)

def handle_br_orientation(v):
    m = re_br_or.findall(v)
    if len(m) > 0:
        f, r = m[0]
        f = int(f)
        r = int(r)
        if f > r:
            return BREAKEND_3P
        elif r > f:
            return BREAKEND_5P
        else:
            return None
    else:
        if v == BREAKEND_3P or v == BREAKEND_5P:
            return v
        else:
            return None

def handle_cg_orientation(or1, or2):
    if or1 == "+":
        o1 = BREAKEND_3P
    else:
        o1 = BREAKEND_5P
    if or2 == "+":  # convert Complete Genomics orientation to standard
        o2 = BREAKEND_5P
    else:
        o2 = BREAKEND_3P
    return o1, o2

def handle_randomize(chrom, pos, randomize, allowed_bp_regions):
    if randomize == False:
        return pos
    if chrom not in allowed_bp_regions:
        return pos
    return random_chr_pos(allowed_bp_regions[chrom])

def canonical_chrom(chrom):
    if chrom.startswith("chr"):
        return chrom[4:]
    elif chrom == "M" or chrom == "25":
        chrom = "MT"
    elif chrom == "23":
        chrom = "X"
    elif chrom == "24":
        chrom = "Y"
    else:
        return chrom

def parse_cg_nodes(f, randomize = False, allowed_bp_regions = None, chroms = None):
    f.readline() # header
    nodes = {}
    next_id = 1
    next_event_id = 1
    for s in f:
        vals = s.strip().split("\t")
        chr1, pos1, or1, chr2, pos2, or2 = vals[:6]
        chr1, chr2 = canonical_chrom(chr1), canonical_chrom(chr2)
        if chroms != None and (chr1 not in chroms or chr2 not in chroms):
            continue
        pos1, pos2 = int(pos1), int(pos2)
        pos1 = handle_randomize(chr1, pos1, randomize, allowed_bp_regions)
        pos2 = handle_randomize(chr2, pos2, randomize, allowed_bp_regions)
        or1, or2 = handle_cg_orientation(or1, or2)
        node1 = Node(next_id, next_event_id, chr1, pos1, or1)
        nodes[next_id] = node1
        next_id += 1
        node2 = Node(next_id, next_event_id, chr2, pos2, or2)
        nodes[next_id] = node2
        next_id += 1
        next_event_id += 1
    return nodes

def parse_breakdancer_nodes(f, randomize = False, allowed_bp_regions = None, chroms = None):
    f.readline() # header
    nodes = {}
    next_id = 1
    next_event_id = 1
    for s in f:
        vals = s.strip().split("\t")
        chr1, pos1, or1, chr2, pos2, or2 = vals[:6]
        chr1, chr2 = canonical_chrom(chr1), canonical_chrom(chr2)
        if chroms != None and (chr1 not in chroms or chr2 not in chroms):
            continue
        pos1 = int(pos1)
        pos2 = int(pos2)
        pos1 = handle_randomize(chr1, pos1, randomize, allowed_bp_regions)
        pos2 = handle_randomize(chr2, pos2, randomize, allowed_bp_regions)
        or1 = handle_br_orientation(or1)
        or2 = handle_br_orientation(or2)
        node1 = Node(next_id, next_event_id, chr1, pos1, or1)
        nodes[next_id] = node1
        next_id += 1
        node2 = Node(next_id, next_event_id, chr2, pos2, or2)
        nodes[next_id] = node2
        next_id += 1
        next_event_id += 1
    return nodes

DELLY_ENDS = {"5to5" : (BREAKEND_5P, BREAKEND_5P),
              "3to3" : (BREAKEND_3P, BREAKEND_3P),
              "3to5" : (BREAKEND_3P, BREAKEND_5P),
              "5to3" : (BREAKEND_5P, BREAKEND_3P)}

DELLY_INV_ENDS = {0 : (BREAKEND_5P, BREAKEND_5P),
                  1 : (BREAKEND_3P, BREAKEND_3P)}

def handle_delly_orientation(svid):
    v = svid.split("_")
    tag = v[1]
    if tag in DELLY_ENDS:
        return DELLY_ENDS[tag]
    else:
        raise Exception("Unknown orientation in DELLY data: %s" % (svid))

def parse_delly_nodes(f, randomize = False, allowed_bp_regions = None, chroms = None):
    nodes = {}
    next_id = 1
    next_event_id = 1
    for s in f:
        if s.startswith("#"):
            continue
        vals = s.strip().split("\t")
        chr1, pos1, chr2, pos2, reads, qual, svid = vals
        chr1, chr2 = canonical_chrom(chr1), canonical_chrom(chr2)
        if chroms != None and (chr1 not in chroms or chr2 not in chroms):
            continue
        pos1, pos2 = int(pos1), int(pos2)
        pos1 = handle_randomize(chr1, pos1, randomize, allowed_bp_regions)
        pos2 = handle_randomize(chr2, pos2, randomize, allowed_bp_regions)
        or1, or2 = handle_delly_orientation(svid)
        node1 = Node(next_id, next_event_id, chr1, pos1, or1)
        nodes[next_id] = node1
        next_id += 1
        node2 = Node(next_id, next_event_id, chr2, pos2, or2)
        nodes[next_id] = node2
        next_id += 1
        next_event_id += 1
    return nodes

def parse_jsva_nodes(f, args, randomize = False, allowed_bp_regions = None, exclude_repeats = set(), chroms = None):
    nodes = {}
    next_id = 1
    next_event_id = 1
    first_bps = {}
    for s in f:
        if s.startswith("#"):
            continue
        bp = unravel.SampleBreakpoint(s)
        if args.omit_multisample and bp.n_somatic != 1:
            continue
        if bp.treads < args.min_tumor_reads:
            continue
        if bp.svid not in first_bps:
            first_bps[bp.svid] = bp
        else:
            bp1 = first_bps[bp.svid]
            bp2 = bp
            if (bp1.repeat[0] in exclude_repeats and bp1.repeat[1] == 0) or (bp2.repeat[0] in exclude_repeats and bp2.repeat[1] == 0):
                continue
            if (bp1.chrom == bp2.chrom and bp1.start > bp2.start) or (bp1.chrom != bp2.chrom and bp1.chrom != bp1.mbp1chrom):
#                sys.stderr.write("Warning: Breakpoint switch %s vs %s\n" % (bp1, bp2)) # this shouldn't happen if unravel.py works as it should
                bp1, bp2 = bp2, bp1
            or1, or2 = JSVA_BREAKEND_ORIENTATION[bp1.or1], JSVA_BREAKEND_ORIENTATION[bp1.or2]
            chr1, chr2 = bp1.mbp1chrom, bp1.mbp2chrom
            chr1, chr2 = canonical_chrom(chr1), canonical_chrom(chr2)
            if chroms != None and (chr1 not in chroms or chr2 not in chroms):
                continue
            pos1, pos2 = int(float(bp1.mbp1pos)), int(float(bp1.mbp2pos))
            pos1 = handle_randomize(chr1, pos1, randomize, allowed_bp_regions)
            pos2 = handle_randomize(chr2, pos2, randomize, allowed_bp_regions)
            n1 = Node(next_id, next_event_id, chr1, pos1, or1)
            n1.reads = (bp1.nreads, bp1.treads)
            n1.gene = bp1.gene
            n1.repeat = bp1.repeat
            nodes[next_id] = n1
            next_id += 1
            n2 = Node(next_id, next_event_id, chr2, pos2, or2)
            n2.reads = (bp2.nreads, bp2.treads) # should be identical to n1.reads
            assert(n1.reads == n2.reads)
            n2.gene = bp2.gene
            n2.repeat = bp2.repeat
            nodes[next_id] = n2
            next_id += 1
            next_event_id += 1
    return nodes

re_delly_tr = re.compile(">Translocation_(.+?)_")
re_delly_dup = re.compile(">Duplication_")
re_delly_inv = re.compile(">Inversion_(\d)_")
re_delly_del = re.compile(">Deletion_")

def handle_delly_filtered_orientation(svid):
    m = re_delly_tr.search(svid)
    if m != None:
        return DELLY_ENDS[m.group(1)]
    if re_delly_dup.search(svid):
        return BREAKEND_5P, "L"
    if re_delly_inv.search(svid):
        return
    if re_delly_del.search(svid):
        return "L", BREAKEND_5P
    raise Exception("Unable to parse orientation: %s\n" % (svid))

def parse_delly_filtered_nodes(f, chroms):
    raise Exception("Filtered DELLY or new DELLY VCF formats are not yet supported")
    nodes = {}
    next_id = 1
    next_event_id = 1
    for s in f:
        if s.startswith("#"):
            continue
        v = s.strip().split("\t")
        svid = v[6]
        if "Translocation" in svid:
            chr1, pos1, chr2, pos2, reads, qual = v[:6]
        else:
            chr1, pos1, pos2, svlen, reads, qual = v[:6]
            chr2 = chr2

        chr1, chr2 = canonical_chrom(chr1), canonical_chrom(chr2)
        if chroms != None and (chr1 not in chroms or chr2 not in chroms):
            continue

        pos1, pos2 = int(pos1), int(pos2)
        or1, or2 = handle_delly_filtered_orientation(svid)
        node1 = Node(next_id, next_event_id, chr1, pos1, or1)
        nodes[next_id] = node1
        next_id += 1
        node2 = Node(next_id, next_event_id, chr2, pos2, or2)
        nodes[next_id] = node2
        next_id += 1
        next_event_id += 1
    return nodes

def build_chr2nodes(nodes):
    chr2nodes = {}
    for nid in nodes:
        n = nodes[nid]
        if n.chrom not in chr2nodes:
            chr2nodes[n.chrom] = []
        chr2nodes[n.chrom].append(n)
    for chrom in chr2nodes:
        chr2nodes[chrom].sort(lambda x, y: cmp(x.pos, y.pos))
    return chr2nodes

def join_nodes(nodes, args):
    """Join sequential breakpoints (nodes) that are closer than args.join_distance to each other. If args.end_matching is set, require that ends match (F-R or R-F). If args.ordered_ends is set, require further that tail precedes head (pos(F)<pos(R))."""
    links = []
    chr2nodes = build_chr2nodes(nodes)
    for chrom in chr2nodes:
        prev_node = None
        for n in chr2nodes[chrom]:
            if prev_node != None:
                diff = n.pos - prev_node.pos
                assert(diff >= 0)
                if n.adj_pos == None or prev_node.adj_pos == None:
                    adj_diff = None
                    proximal = diff <= args.join_distance
                else:
                    adj_diff = n.adj_pos - prev_node.adj_pos
                    proximal = adj_diff <= args.join_distance

                order = args.ordered_ends == False or prev_node.end == BREAKEND_3P and n.end == BREAKEND_5P
                match = args.end_matching == False or n.end != prev_node.end

                if proximal and match and order:
#                    print "JOIN", prev_node, n, diff, adj_diff
                    pass

                    links.append((prev_node.node_id, n.node_id, diff, adj_diff))
            prev_node = n
    return links

def parse_breakdancer(f, cnv, randomize, allowed_bp_regions, args):
    """Parse breakpoints in BreakDancer format. Orientations may be in original (X+Y-) or called (F/R) format."""
    nodes = parse_breakdancer_nodes(f, randomize, allowed_bp_regions, args.chroms)
    contract(nodes, cnv, args)
    links = join_nodes(nodes, args)
    return nodes, links

def parse_cg(f, cnv, randomize, allowed_bp_regions, args):
    nodes = parse_cg_nodes(f, randomize, allowed_bp_regions, args.chroms)
    contract(nodes, cnv, args)
    links = join_nodes(nodes, args)
    return nodes, links

def parse_delly(f, cnv, randomize, allowed_bp_regions, args):
    nodes = parse_delly_nodes(f, randomize, allowed_bp_regions, args.chroms)
    contract(nodes, cnv, args)
    links = join_nodes(nodes, args)
    return nodes, links

def parse_delly_filtered(f, cnv, randomize, allowed_bp_regions, args):
    nodes = parse_delly_filtered_nodes(f, randomize, allowed_bp_regions, args.chroms)
    contract(nodes, cnv, args)
    links = join_nodes(nodes, args)
    return nodes, links

def parse_annotated(f, randomize, allowed_bp_regions, args):
    next_id = 1
    expecting_head = False
    prev_node = None
    nodes = {}
    end_links = []
    for s in f:
        vals = s.strip().split("\t")
        if len(vals) < 4:
            continue
        event, chrom, pos, end = vals[:4]
        chrom = canonical_chrom(chrom)
        if args.chroms != None and chrom not in args.chroms:
            continue
        event = int(event)
        pos = handle_randomize(chrom, int(pos), randomize, allowed_bp_regions)

        node = Node(next_id, event, chrom, pos, end)
        nodes[next_id] = node
        next_id += 1

        if end == "H" and expecting_head:
            # link to previous tail end
            if prev_node.chrom == node.chrom:
                diff = node.pos - prev_node.pos
            else:
                diff = "NA"
            end_links.append((prev_node.node_id, node.node_id, diff))

        if end == "T":
            prev_node = node
            expecting_head = True

    return nodes, end_links

def parse_jsva(f, cnv, randomize, allowed_bp_regions, args):
    """Parse joint structural variant analysis result file format generated by unravel.py.
    Assume that the input contains breakpoints from only a single tumor."""
    nodes = parse_jsva_nodes(f, args, randomize, allowed_bp_regions, args.exclude_repeat, args.chroms)
    contract(nodes, cnv, args)
    links = join_nodes(nodes, args)
    return nodes, links

def add_link(G, u, v, w, diff, adj_diff):
    if u not in G:
        G[u] = {}
    G[u][v] = (w, diff, adj_diff)
    if v not in G:
        G[v] = {}
    G[v][u] = (w, diff, adj_diff)

def build_graph(nodes, end_links):
    # graph is constructed by linking breakpoints that
    # a) belong to the same head-tail pair; or
    # b) are associated with the same event

    G = {}
    for u, v, diff, adj_diff in end_links:
        add_link(G, u, v, LINK_END, diff, adj_diff)

    events = {}
    for nid in nodes:
        n = nodes[nid]
        if n.event_id not in events:
            events[n.event_id] = [n.node_id]
        else:
            events[n.event_id].append(n.node_id)

    for e in events:
        for u in events[e]:
            for v in events[e]:
                if u != v:
                    add_link(G, u, v, LINK_EVENT, None, None)

    return G

def fnum(x):
    return locale.format("%d", x, grouping = True)

def graph_to_dot(G, nodes, comp, of, clusters = True, chr_colors = None, min_n = 0):
    of.write("graph rearrangements {\n")
    cl2node = {}
    for u in comp:
        c = comp[u]
        if c not in cl2node:
            cl2node[c] = set()
        cl2node[c].add(u)

    n_nodes = 0
    for c in cl2node:
        if len(cl2node[c]) < min_n:
            continue
        n_nodes += len(cl2node[c])
        if clusters:
            of.write("subgraph cluster_%s {\n" % (c))
        for u in cl2node[c]:
            n = nodes[u]
            if n.chrom in chr_colors:
                bgcolor, fgcolor = chr_colors[n.chrom]
            else:
                bgcolor, fgcolor = "white", "black"
            if fgcolor == bgcolor:
                bgcolor, fgcolor = "white", "black"
            pos_str = locale.format("%d", n.pos, grouping = True)
            gene_str = repeat_str = reads_str = ""
            if n.gene[0] != None:
                gene_str = "%s (%s)" % (n.gene[0], fnum(n.gene[1]))
            if n.repeat[0] != None:
                repeat_str = "%s (%s)" % (n.repeat[0], fnum(n.repeat[1]))
            if n.reads[0] != None:
                reads_str = "NR=%s TR=%s" % (n.reads[0], n.reads[1])

            label = "%s event=%s chr%s\\npos=%s end=%s\\n%s\\n%s\\n%s" % (n.node_id, n.event_id, n.chrom, pos_str, n.end, gene_str, repeat_str, reads_str)
            of.write("  %s [label=\"%s\",style=\"filled\",fillcolor=\"%s\",fontcolor=\"%s\"];\n" % (n.node_id, label, bgcolor, fgcolor))
            for v in G[u]:
                if u <= v:
                    w, diff, adj_diff = G[u][v]
                    if w == LINK_EVENT:
                        col = "blue"
                    elif w == LINK_END:
                        col = "red"
                    else:
                        assert(0) # invalid link type
                    if adj_diff != None:
                        label = "%s->%s" % (fnum(diff), fnum(adj_diff))
                    elif diff != None:
                        label = "%s" % (fnum(diff))
                    else:
                        label = ""
                    of.write("   %s -- %s [color=\"%s\",label=\" %s\"];\n" % (u, v, col, label))
        if clusters:
            of.write("}\n")

    of.write("}\n")
#    print "%d nodes total" % (n_nodes)

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

def yesno(b):
    if b == True:
        return "Yes"
    elif b == False:
        return "No"
    else:
        assert(0)

def write_results(o, nodes, end_links, comp, cyclic_comp, complex_calls):
    o.write("#Breakpoint event cluster associations\n")
    o.write("#Breakpoint Event Chr Pos End EventCluster IsEventCyclic IsEventComplex\n")
    for nid in nodes:
        n = nodes[nid]
        c = comp[nid]
        cyclic = yesno(cyclic_comp[c])
        chrt = yesno(complex_calls[c])
        o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (n.node_id, n.event_id, n.chrom, n.pos, n.end, c, cyclic, chrt))

def read_cnv_data_varscan(f, args):
    cnv = {}
    for s in f:
        if s.startswith("#"):
            continue
        chrom, start, end, numsegs, logcnratio, call = s.strip().split("\t")
        chrom = canonical_chrom(chrom)
        start, end = int(start), int(end)
        if chrom not in cnv:
            cnv[chrom] = []
        cnv[chrom].append((start, end, 1, float(logcnratio)))
    for chrom in cnv:
        cnv[chrom].sort(lambda x, y: cmp(x[0], y[0]))
    return cnv

def read_cnv_data_raw(f, args):
    cnv = {}
    change_sum = 0.0
    seg_count = 0
    for s in f:
        vals = s.strip().split("\t")
        chrom, start, end, support, change = vals

        chrom = canonical_chrom(chrom)
        start = int(start)
        end = int(end)
        support = int(support)
        change = float(change)
        if chrom not in cnv:
            cnv[chrom] = []
        cnv[chrom].append((start, end, support, change))
        seg_count += 1
        change_sum += change

    # normalize by subtracting mean log ratio
    if not args.skip_normalize_cnv:
        mean_change = change_sum / seg_count
        for chrom in cnv:
            norm_cnv = []
            for start, end, support, change in cnv[chrom]:
                norm_cnv.append((start, end, support, change - mean_change))
            cnv[chrom] = norm_cnv

    # ensure that breakpoints are ordered by start position
    for chrom in cnv:
        cnv[chrom].sort(lambda x, y: cmp(x[0], y[0]))

    return cnv

def read_cnv_data_called(f, args):
    cnv = {}
    f.readline()
    for s in f:
        chrom, start, end, sid, cncall, size = s.strip().split("\t")
        chrom = canonical_chrom(chrom)
        start = int(start)
        end = int(end)
        if chrom not in cnv:
            cnv[chrom] = []
        if cncall == "Deletion" or cncall == "Double Deletion":
            change = args.deletion
        else:
            change = 0
        support = 1
        cnv[chrom].append((start, end, support, change))
    for chrom in cnv:
        cnv[chrom].sort(lambda x, y: cmp(x[0], y[0]))
    return cnv

def read_cnv_data(f, args):
    if args.cnv_data_format == CNV_DATA_FORMAT_RAW:
        return read_cnv_data_raw(f, args)
    elif args.cnv_data_format == CNV_DATA_FORMAT_CALLED:
        return read_cnv_data_called(f, args)
    elif args.cnv_data_format == CNV_DATA_FORMAT_VARSCAN:
        return read_cnv_data_varscan(f, args)
    else:
        assert(0) # invalid cnv data format

def contract(nodes, cnv, args):
    """Adjust breakpoint positions by considering deletions of chromosomal regions. Deletions are inferred from CNV data."""
    if cnv == None:
        return

    # We set up a list of breakpoints of two types:
    # 1) breakpoints from paired-end data (nodes)
    chr2node = {}
    for nid in nodes:
        n = nodes[nid]
        if n.chrom not in chr2node:
            chr2node[n.chrom] = []
        cn = ContractNode()
        cn.pos = n.pos
        cn.bptype = ContractNode.NODE_BP
        cn.node_id = n.node_id
        chr2node[n.chrom].append(cn)

    # 2) breakpoints from CNV data indicating deletions (cnv)
    for chrom in cnv:
        for (start, end, support, change) in cnv[chrom]:
            if change >= args.deletion:
                # call non-deletion, ignore
                continue
            cn1 = ContractNode()
            cn1.pos = start
            cn1.bptype = ContractNode.NODE_CNV_START
            cn2 = ContractNode()
            cn2.pos = end
            cn2.bptype = ContractNode.NODE_CNV_END
            if chrom not in chr2node:
                chr2node[chrom] = []
            chr2node[chrom].append(cn1)
            chr2node[chrom].append(cn2)

    for chrom in chr2node:
        chr2node[chrom].sort(lambda x, y: cmp(x.pos, y.pos))
        del_region_start = None
        pos_adj = 0
        for cn in chr2node[chrom]:
            if cn.bptype == ContractNode.NODE_CNV_START:
                del_region_start = cn.pos
                to_be_deleted = []
            elif cn.bptype == ContractNode.NODE_CNV_END:
                pos_adj -= cn.pos - del_region_start
                del_region_start = None
                for cnd in to_be_deleted:
                    if cn.pos - cnd.pos > args.deletion_margin:
                        nodes[cnd.node_id].adj_pos = None  # delete only candidate node is distant enough from cnv segment start and end
#                       print "Deleted breakpoint - inside deletion region"

            elif cn.bptype == ContractNode.NODE_BP:
                if del_region_start != None:
                    # inside deletion region
                    if args.deletion_strategy == STRATEGY_DEL_BP_DELETE:
                        # delete bp
                        if cn.pos - del_region_start > args.deletion_margin:
                            to_be_deleted.append(cn)  # delay the deletion until we can check distance from cnv segment end
#                            nodes[cn.node_id].adj_pos = None
                    elif args.deletion_strategy == STRATEGY_DEL_BP_REPOSITION:
                        nodes[cn.node_id].adj_pos = del_region_start + pos_adj
#                        print "Repositioned breakpoint %s -> %s - inside deletion region" % (nodes[cn.node_id].pos, nodes[cn.node_id].adj_pos)
                else:
                    nodes[cn.node_id].adj_pos = nodes[cn.node_id].pos + pos_adj
#                    print "Adjusted", nodes[cn.node_id].pos, "->", nodes[cn.node_id].adj_pos
            else:
                assert(0) # invalid cn.bptype

    return nodes

def adjusted_end_links(end_links, nodes):
    E = []
    for u, v, diff in end_links:
        if nodes[v].adj_pos != None and nodes[u].adj_pos != None:
            ad = nodes[v].adj_pos - nodes[u].adj_pos
        else:
            ad = None
        E.append((u, v, diff, ad))
    return E

def write_end_links(o, nodes, end_links):
    o.write("#Breakpoints linked by proximity\n")
    o.write("#Distance: distance between two breakpoints\n")
    o.write("#AdjustedDistance: breakpoint distance after calling deletions in copy number data (None if no CNV data used)\n")
    o.write("#Breakpoint1 Chr1 Pos1 End1 Breakpoint2 Chr2 Pos2 End2 Distance AdjustedDistance\n")
    for u, v, diff, adj_diff in end_links:
        prev_node = nodes[u]
        node = nodes[v]
        o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (prev_node.node_id, prev_node.chrom, prev_node.pos, prev_node.end, node.node_id, node.chrom, node.pos, node.end, diff, adj_diff))

def build(f, cnv, args, randomize = False, allowed_bp_regions = None):
    if args.input_type == INPUT_ANNOTATED:
        nodes, end_links = parse_annotated(f, randomize, allowed_bp_regions, args)
        contract(nodes, cnv, args)
        end_links = adjusted_end_links(end_links, nodes)
        return nodes, end_links
    elif args.input_type in PARSERS:
        return PARSERS[args.input_type](f, cnv, randomize, allowed_bp_regions, args)
    else:
        print "Unsupported input type %s" % (args.input_type)
        sys.exit(1)

def cluster_sort_fun(x, y):
    if x[0] != y[0]:
        return cmp(y[0], x[0])
    else:
        return cmp(y[2], x[2])

def write_summary(o, nodes, comp, cyclic_comp, complex_calls, p, eq, ts, args):
    c2n = {}
    for u in comp:
        c = comp[u]
        if c not in c2n:
            c2n[c] = set()
        c2n[c].add(u)
    items = []
    nc = 0
    for c in c2n:
        chroms = {}
        total = 0
        for u in c2n[c]:
            n = nodes[u]
            if n.chrom not in chroms:
                chroms[n.chrom] = 1
            else:
                chroms[n.chrom] += 1
            total += 1
        chrs = chroms.keys()
        chrs.sort()
        intra = yesno(len(chrs) == 1)
        if complex_calls[c]:
            complex = "Yes"
            nc += 1
        else:
            complex = "No"
        cyclic = yesno(cyclic_comp[c])
        items.append((len(chroms), ",".join(map(str, chrs)), total, cyclic, complex, intra))

    items.sort(cluster_sort_fun)
    comp_ix = 1
    o.write("#Number of complex events: %d\n" % (nc))
    o.write("#Complex event test p %s %s (%d rounds)\n" % (eq, p, args.rounds))
    o.write("#EventCluster NumChrs Chrs NumBreakpoints Intrachromosomal CyclicGraph ComplexCall\n")
    for nchroms, chrs, total, cyclic, complex, intra in items:
        o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (comp_ix, nchroms, chrs, total, intra, cyclic, complex))
        comp_ix += 1

def convert_color(s):
    r, g, b = map(int, s.strip().split(","))
    if r + g + b < 205:
        fgcol = "white"
    else:
        fgcol = "black"
    return ("#%.2x%.2x%.2x" % (r, g, b), fgcol)

def load_colors(f):
    colors = {}
    for s in f:
        chrom, col = s.strip().split("\t")
        colors[chrom.strip()] = convert_color(col)
    return colors

def component_sizes(comp):
    comp_size = {}
    for u in comp:
        c = comp[u]
        if c not in comp_size:
            comp_size[c] = 1
        else:
            comp_size[c] += 1
    return comp_size

def permute_breakpoint_pairs(f):
    # assume Breakdancer raw format
    m1 = []
    m2 = []
    f.readline()
    for s in f:
        chr1, pos1, or1, chr2, pos2, or2 = s.strip().split("\t")
        m1.append((chr1, pos1, or1))
        m2.append((chr2, pos2, or2))
    R = ""
    random.shuffle(m2)
    i = 0
    for chr1, pos1, or1 in m1:
        chr2, pos2, or2 = m2[i]
        i += 1
        R += "%s\t%s\t%s\t%s\t%s\t%s\n" % (chr1, pos1, or1, chr2, pos2, or2)
    return R

def update_region(chr1, p, chr_regions):
    if chr1 not in chr_regions:
        chr_regions[chr1] = [p, p]
    else:
        chr_regions[chr1][0] = min(chr_regions[chr1][0], p)
        chr_regions[chr1][1] = max(chr_regions[chr1][1], p)

def random_chr_pos(chr_regions):
    region_size, regions = chr_regions
    r = random.randint(0, region_size)
    pos = 0
    for start, end in regions:
        if end - start + pos > r:
            pos = start + r - pos
            break
        else:
            pos += end - start
    return pos

def permute_breakpoints(f, allowed_bp_regions):
    # assume Breakdancer raw or Complete Genomics format
    M = []
    f.readline()
    for s in f:
        vals = s.strip().split("\t")
        chr1, pos1, or1, chr2, pos2, or2 = vals[:6]
        pos1 = int(pos1)
        pos2 = int(pos2)
        M.append(vals)
    R = ""
    for vals in M:
        chr1, pos1, or1, chr2, pos2, or2 = vals[:6]
        p1 = random_chr_pos(allowed_bp_regions[chr1])
        p2 = random_chr_pos(allowed_bp_regions[chr2])
        R += "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr1, p1, or1, chr2, p2, or2, "\t".join(vals[6:]))
    return R

def call_complex(comp_sizes, cyclic_comp, args):
    t = 0
    calls = {}
    for c in comp_sizes:
        cn = comp_sizes[c]
        call = cn >= args.min_complex_size and not (cn == 4 and cyclic_comp[c] == True)
        calls[c] = call
    return calls

def test_statistic(comp_sizes, complex_calls, args):
    """Count the number of event clusters that are called complex:
    a complex event cluster has at least n=min_complex_size breakpoints and,
    in case n >= 4, size 4 cluster has to be non-cyclic (otherwise balanced translocations
    would be accounted as complex events)."""
    if args.statistic == STATISTIC_EVENT_COUNT:
        return len(filter(lambda x: x, complex_calls.values())) # number of complex events
    elif args.statistic == STATISTIC_EVENT_BREAKPOINT_COUNT:
        bp = 0
        for c in complex_calls:
            if complex_calls[c]:
                bp += comp_sizes[c]
        return bp
    else:
        raise Exception("Unknown test statistic %d" % (args.statistic))

def visualize_event_clusters(G, nodes, comp, fname, chr_colors, args):
    locale.setlocale(locale.LC_ALL, "en_US.UTF-8")
    vf = tempfile.NamedTemporaryFile()
    vf2 = tempfile.NamedTemporaryFile()
    graph_to_dot(G, nodes, comp, vf, clusters = False, chr_colors = chr_colors, min_n = args.min_cluster_size)
    graph_to_dot(G, nodes, comp, vf2, clusters = True, chr_colors = chr_colors, min_n = args.min_cluster_size)
    vf.flush()
    vf2.flush()
#    subprocess.call("dot %s -o %s.png -Tpng" % (vf.name, fname), shell = True)
    subprocess.call("dot %s -o %s.clusters.png -Tpng" % (vf2.name, fname), shell = True)
    vf.close()
    vf2.close()

def randomize(f, cnv, allowed_bp_regions, chr_colors, args):
#    rdata = StringIO.StringIO(permute_breakpoints(f, allowed_bp_regions))
    nodes, end_links = build(f, cnv, randomize = True, allowed_bp_regions = allowed_bp_regions, args = args)
    G = build_graph(nodes, end_links)
    comp = find_components(G)
    comp_sizes = component_sizes(comp)
    cyclic_comp = is_cyclic(G, comp)

    if args.visualize_permutations:
        visualize_event_clusters(G, nodes, comp, "%srandom%d" % (args.output, random.randint(1, 100000)), chr_colors, args)

    return test_statistic(comp_sizes, call_complex(comp_sizes, cyclic_comp, args), args)

def read_bed(f):
    regions = {}
    for s in f:
        if s.startswith("#"):
            continue
        vals = s.strip().split("\t")
        chrom, start, end = vals[:3]
        start, end = int(start), int(end)
        if chrom not in regions:
            regions[chrom] = []
        regions[chrom].append((start, end))
    return regions

def min_max_regions(regions):
    min_max = {}
    for chrom in regions:
        mn = sys.maxint
        mx = -sys.maxint
        for start, end in regions[chrom]:
            if start < mn:
                mn = start
            if end > mx:
                mx = end
        min_max[chrom] = (mn, mx)
    return min_max

def combine_regions(min_max, excluded_regions = None):
    regions = {}
    for chrom in min_max:
        mn, mx = min_max[chrom]
        assert(mx >= mx)

        regions[chrom] = [0, []]  # chr -> size, list of (start, end)
        begin = mn
        n = 0
        if excluded_regions != None and chrom in excluded_regions:
            for (start, end) in excluded_regions[chrom]:
                regions[chrom][1].append((begin, start - 1))
                n += start - begin - 1
                begin = end + 1
            regions[chrom][1].append((begin, mx))
            regions[chrom][0] = n
            assert(n >= 0)
        else:
            regions[chrom] = [mx - mn, [(mn, mx)]]

    return regions

def plot_permutation_test(ts, rts, rounds, title, outfn):
    plt.figure()
    plt.hist(rts, bins = max(20, int(1.0 * rounds / 10)), color = "blue")
    plt.axvline(ts, color = "red")
    plt.title(title)
    plt.savefig(outfn, bbox_inches = "tight", dpi = 300)

def parse_chr_len(f):
    l = {}
    for s in f:
        chrom, slen = s.strip().split("\t")
        l[chrom] = (0, int(slen))
    return l

def perform_statistics(data, cnv, allowed_bp_regions, chr_colors, complex_calls, comp_sizes, args):
    of = open("%s.test" % (args.output), "w")
    of.write("#Sample n_complex n_rounds p_eq p_value rnd_50%% rnd_5%% rnd_95%%\n")


    ts = test_statistic(comp_sizes, complex_calls, args)

    if ts == 0:
        of.write("%s\n" % ("\t".join(map(str, [args.output, ts, args.rounds, "=", 1.0, np.nan, np.nan, np.nan]))))
        return 1.0, "=", 0

    if args.rounds == 0:
        of.write("%s\n" % ("\t".join(map(str, [args.output, ts, args.rounds, "=", "N/A", np.nan, np.nan, np.nan]))))
        return "N/A", "=", ts

    # randomize
    rts = []
    for i in range(args.rounds):
        sys.stderr.write(".")
        sys.stderr.flush()
        rts.append(randomize(StringIO.StringIO(data), cnv, allowed_bp_regions, chr_colors, args))
    rts.sort()
    sys.stderr.write("\n")

    if args.null:
        on = open("%s.null" % (args.output), "w")
        on.write("#Output of \"%s\" on %s\n" % (" ".join(sys.argv), datetime.datetime.now()))
        on.write("#Complex event test, null distribution\n")
        for x in rts:
            on.write("%s\n" % (x))
        on.close()

    l = bisect.bisect_left(rts, ts)
    if len(rts) == l:
        p = 1.0 / l
        eq = "<"
    else:
        p = (1.0 * len(rts) - l) / len(rts)
        eq = "="

    of.write("%s\n" % ("\t".join(map(str, [args.output, ts, args.rounds, eq, p, rts[int(0.5 * len(rts))], rts[int(0.05 * len(rts))], rts[int(0.95 * len(rts))]]))))

    plot_permutation_test(ts, rts, args.rounds, "%s (n=%d, p%s%f)" % (args.output, ts, eq, p), "%s.test.png" % (args.output))

    return p, eq, ts

def __main(args):
    if args.input != None and args.input != "-":
        if args.input.endswith(".gz"):
            f = gzip.open(args.input)
        else:
            f = open(args.input)
    else:
        f = sys.stdin
        print "Reading input from stdin"
    if args.output != None and args.output != "-":
        o = open("%s.bp" % (args.output), "w")
        odiff = open("%s.proximal" % (args.output), "w")
        osum = open("%s.txt" % (args.output), "w")
        s = "#Output of \"%s\" on %s\n" % (" ".join(sys.argv), datetime.datetime.now())
        o.write(s)
        odiff.write(s)
        osum.write(s)
    else:
        o = odiff = osum = sys.stdout
    if args.cnv != None:
        if os.path.exists(args.cnv):
            if args.cnv.endswith(".gz"):
                cnv = read_cnv_data(gzip.open(args.cnv), args)
            else:
                cnv = read_cnv_data(open(args.cnv), args)
        else:
            cnv = None
            print "Copy number data file %s does not exist - skipping" % (args.cnv)
    else:
        cnv = None
    if args.disable_colors == False:
        chr_colors = load_colors(open(args.chrom_colors))
    else:
        chr_colors = {}

    chr_regions = parse_chr_len(open(args.ref_chr_len))

    if args.exclude_region != None:
        chr_exclude = read_bed(open(args.exclude_region))
    else:
        chr_exclude = None
    allowed_bp_regions = combine_regions(chr_regions, chr_exclude)

    data = f.read()
    bpdata = StringIO.StringIO(data)

    nodes, end_links = build(bpdata, cnv, args)
    G = build_graph(nodes, end_links)
    comp = find_components(G)
    comp_sizes = component_sizes(comp)
    cyclic_comp = is_cyclic(G, comp)

    sizes = comp_sizes.values()
    sizes.sort(lambda x, y: cmp(y, x))
#    print "Component size distribution:", sizes
#    print "Cyclic components:", cyclic_comp.values()
    calls = call_complex(comp_sizes, cyclic_comp, args)
    #print "Complex calls:", calls.values()

    sys.stdout.write("Complex event status: ")
    sys.stdout.flush()
    p, eq, ts = perform_statistics(data, cnv, allowed_bp_regions, chr_colors, calls, comp_sizes, args)
    print "statistic = %d  (p %s %s, %d rounds)" % (ts, eq, p, args.rounds)

    print "===\t%s\t%s\t%s\t%s" % (args.input, ts, eq, p)

    sys.stdout.write("Writing results... ")
    sys.stdout.flush()
    write_end_links(odiff, nodes, end_links)
    write_results(o, nodes, end_links, comp, cyclic_comp, calls)
    write_summary(osum, nodes, comp, cyclic_comp, calls, p, eq, ts, args)
    sys.stdout.write("done\n")

    if args.visualize:
        if args.output != None and args.output != "-":
            fn = args.output
        else:
            fn = "events"
        sys.stdout.write("Plotting event clusters... ")
        sys.stdout.flush()
        visualize_event_clusters(G, nodes, comp, fn, chr_colors, args)
        sys.stdout.write("done\n")

PARSERS = {INPUT_BREAKDANCER_FR : parse_breakdancer,
           INPUT_BREAKDANCER_RAW : parse_breakdancer,
           INPUT_CG_RAW : parse_cg,
           INPUT_DELLY : parse_delly,
           INPUT_DELLY_FILTERED : parse_delly_filtered,
           INPUT_JSVA : parse_jsva}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument("-i", "--input", help = "Input breakpoints", metavar = "FN")
    parser.add_argument("-c", "--cnv", help = "Input copy number data (optional)", metavar = "FN")
    parser.add_argument("-o", "--output", help = "Output filename base", metavar = "FN")
    parser.add_argument("-z", "--visualize", help = "Write a visualization of the event graph to file", action = "store_true")
    parser.add_argument("-t", "--input-type", help = "Breakpoint list type (%s=annotated, %s=BreakDancer FR, %s=BreakDancer raw, %s=CG, %s=DELLY raw, %s=DELLY filtered, %s=JSVA/unravel; default=%s)" % (INPUT_ANNOTATED, INPUT_BREAKDANCER_FR, INPUT_BREAKDANCER_RAW, INPUT_CG_RAW, INPUT_DELLY, INPUT_DELLY_FILTERED, INPUT_JSVA, DEFAULT_INPUT_TYPE), metavar = "int", type = int, required = True)
    parser.add_argument("-j", "--join-distance", help = "Maximum distance between linked breakpoints (default=%d)" % (DEFAULT_JOIN_DISTANCE), type = int, default = DEFAULT_JOIN_DISTANCE, metavar = "int")
    parser.add_argument("-d", "--deletion", help = "CNV deletion threshold (default=%s)" % (DEFAULT_CNV_DEL_THRESHOLD), type=float, default=DEFAULT_CNV_DEL_THRESHOLD, metavar = "float")
    parser.add_argument("--skip-normalize-cnv", help = "Do not normalize CNV data", default = False, action = "store_true")
    parser.add_argument("--deletion-strategy", help = "Breakpoints within deleted regions: 1=reposition, 2=delete (default=%s)" % (DEFAULT_DEL_BP_STRATEGY), default = DEFAULT_DEL_BP_STRATEGY, type = int, metavar = "int")
    parser.add_argument("--deletion-margin", help = "If using --deletion-strategy=2, specify safety margin in bp [%(default)d]", default = 10000, type = int, metavar = "int")
    parser.add_argument("--end-matching", help = "Require that joined breakpoints have matching ends", action = "store_true", default = False)
    parser.add_argument("--ordered-ends", help = "Require that tail end position precedes head end position when joining breakpoints", action = "store_true", default = False)
    parser.add_argument("--disable-colors", help = "Do not color chromosomes according to UCSC coloring scheme", action = "store_true", default = False)
    parser.add_argument("-m", "--min-cluster-size", help = "Minimum event cluster size to visualize [%(default)d]", default = 4, type = int)
    parser.add_argument("--min-complex-size", help = "Minimum event cluster size to call complex cluster (default=%d)" % (DEFAULT_MIN_COMPLEX_SIZE), default = DEFAULT_MIN_COMPLEX_SIZE, type = int)
    parser.add_argument("-r", "--rounds", help = "Randomization rounds (default=%d)" % (DEFAULT_ROUNDS), type = int, default = DEFAULT_ROUNDS)
    parser.add_argument("--cnv-data-format", help = "Copy number data format: %d=chr/start/end/len/change, %d=called, %d=VarScan [%d]" % (CNV_DATA_FORMAT_RAW, CNV_DATA_FORMAT_CALLED, CNV_DATA_FORMAT_VARSCAN, CNV_DATA_FORMAT_VARSCAN), default = CNV_DATA_FORMAT_VARSCAN, type = int)
    parser.add_argument("-x", "--exclude-region", help = "Exclude genomic regions (BED file)", metavar = "BED")
    parser.add_argument("-p", "--exclude-repeat", help = "Exclude RepeatMasker region [%(default)s]", default = ",".join(DEFAULT_EXCLUDE_REPEATS))
    parser.add_argument("--chroms", help = "Chromosomes to consider, comma-separated, \"-\" allows all [1..22,X,Y]", default = DEFAULT_CHROMS)
    parser.add_argument("-s", "--statistic", help = "Complex event test statistic: %d=event count, %d=event breakpoint count (default=%d)" % (STATISTIC_EVENT_COUNT, STATISTIC_EVENT_BREAKPOINT_COUNT, DEFAULT_STATISTIC), type = int, default = DEFAULT_STATISTIC)
    parser.add_argument("-n", "--null", help = "Write null distribution", action = "store_true")
    parser.add_argument("--omit-multisample", help = "Only use events appearing in a single tumor (JSVA input only)", action = "store_true")
    parser.add_argument("--min-tumor-reads", help = "Minimum tumor read support (JSVA input only) [%(default)d]", type = int, default = 0)
    parser.add_argument("--visualize-permutations", help = "Draw permuted event clusters", action = "store_true")
    parser.add_argument("--chrom-colors", default="annotation/UCSC_chr_colors.txt")
    parser.add_argument("--ref-chr-len", default="reference/hs37d5_viral.fa.seqlen.bedtools")

    args = parser.parse_args()
    if args.chroms == "-":
        args.chroms = None
    elif type(args.chroms) == str:
        args.chroms = args.chroms.split(",")
    __main(args)
