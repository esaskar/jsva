#!/usr/bin/env python

import sys, argparse, gzip, errno

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input")
    p.add_argument("-o", "--output")
    args = p.parse_args()

    if args.input == None:
        f = sys.stdin
    else:
        if args.input.endswith(".gz"):
            f = gzip.open(args.input)
        else:
            f = open(args.input)
    if args.output == None:
        o = sys.stdout
    else:
        if args.output.endswith(".gz"):
            o = gzip.open(args.output, "w")
        else:
            o = open(args.output, "w")

    try:
        for s in f:
            if s.startswith("#"):
                continue
            elif s.startswith("Chrom"):
                o.write(s)
                continue
            line = s.strip().split("\t")
            if line[5] == "tr":
                chrom1, chrom2 = line[0].split("/")
            else:
                chrom1 = chrom2 = line[0]
            o.write("%s\t%s\t%s\t%s\n" % (chrom1, line[1], line[1], "\t".join(line[3:])))
            o.write("%s\t%s\t%s\t%s\n" % (chrom2, line[2], line[2], "\t".join(line[3:])))
#            svid += 1
    except IOError, e:
        if e.errno != errno.EPIPE:
            raise
