#!/usr/bin/env python

import sys, argparse, os, gzip, errno
import circos
from breakpoint import Breakpoint, BreakpointSet
from unravel import SampleBreakpoint

SVTYPE_COLOR = {"del" : "red", "inv" : "blue", "dup" : "green", "tr" : "orange"}

SOMATIC_COLORS = [(0.9, "red"), 
                  (0.05, "orange"),
                  (0.04, "green"),
                  (0.03, "blue"),
                  (0.02, "purple"),
                  (0.01, "grey")]

def get_link_color_thickness(bp, total_somatic):
    if bp.ngermline > 0:
        return "vlgrey", 1
    else:
        for threshold, col in SOMATIC_COLORS:
            f = 1.0 * bp.nsomatic / total_somatic
            if f > threshold:
#                print bp.nsomatic, f, col, min(10, int(f * 10))
                return col, max(1, min(10, int(f * 10)))
        return "vlgrey", 1

def get_link_z(bp):
    z = 1000 - bp.nsomatic
    return z

def main(args):
    if args.input == None:
        inf = sys.stdin
    else:
        try:
            inf = gzip.open(args.input)
        except:
            inf = open(args.input)
    if args.header == None:
        if args.input != None:
            if args.input.endswith(".tsv.gz"):
                n = len(".tsv.gz")
            else:
                n = len(".tsv")
            args.header = "%s.header.txt" % (args.input[:-n])
        else:
            sys.stderr.write("Specify header name when using stdin\n")
            sys.exit(2)

    bps = BreakpointSet(inf, args.header)

    # Plot tumor type specific circos plots
    links = []
    histogram = []
    text = []
    genes = {}
    c = 0
    for bp in bps:
        if bp.svtype == Breakpoint.SV_TR or bp.length > args.min_link_length:
            #            links.append((chr1, pos1, chr2, pos2, color, thickness, label))
            if bp.ngermline > 0:
                continue
            if bp.nsomatic < 3:
                continue
            chr1 = bp.data[bps.get_column(BreakpointSet.IX_BP1_CHR)]
            pos1 = bp.data[bps.get_column(BreakpointSet.IX_BP1_POS)]
            chr2 = bp.data[bps.get_column(BreakpointSet.IX_BP2_CHR)]
            pos2 = bp.data[bps.get_column(BreakpointSet.IX_BP2_POS)]
            color, thickness = get_link_color_thickness(bp, bps.num_somatic)
            z = get_link_z(bp)
            links.append(circos.Link(chr1, int(float(pos1)), chr2, int(float(pos2)), color = color, thickness = thickness, z = z)
#            c += 1
#            if c > 1000:
#                break
    
    circos.plot("%s/WGS" % (args.output), links, histogram, text)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input")
    p.add_argument("-e", "--header")
    p.add_argument("-o", "--output", required = True)
    p.add_argument("--min-link-length", default = 10000, type = int, help = "[%(default)s]")
    p.add_argument("--gene-dist", default = 10000, type = int, help = "[%(default)s]")
    p.add_argument("--min-reads-label", default = 5, type = int, help = "Minimum number of reads to annotate SV with gene name [%(default)s]")
    args = p.parse_args()
    main(args)
