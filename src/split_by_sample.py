#!/usr/bin/env python2
"""
Split a file from unravel.py by sample.
"""

import sys, argparse, gzip
import unravel

def __main(args):
    if args.input == None or args.input == "-":
        sys.stderr.write("Reading from stdin...\n")
        f = sys.stdin
    elif args.input.endswith(".gz"):
        f = gzip.open(args.input)
    else:
        f = open(args.input)

    sample_f = {}
        
    for s in f:
        if s.startswith("#"):
            continue
        if s.startswith("Chrom"):
            continue
        bp = unravel.SampleBreakpoint(s)
        if bp.sample not in sample_f:
            sample_f[bp.sample] = gzip.open("%s/%s.gz" % (args.outdir, bp.sample), "wb")
#            sys.stderr.write("%s\n" % (bp.sample))
        sample_f[bp.sample].write(s)
        
    for s in sample_f:
        sample_f[s].close()

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input")
    p.add_argument("-o", "--outdir", required = True)
    args = p.parse_args()
    __main(args)