#!/usr/bin/env python

import sys

for s in sys.stdin:
    if s.startswith("#"):
        sys.stdout.write(s)
    else:
        line = s.strip().split("\t")
        if line[4] == "tr":
            chrom1, chrom2 = line[0].split("/")
            sys.stdout.write("%s\t%s\t%s\t%s\n" % (chrom1, line[1], line[1], "\t".join(line[3:])))
            sys.stdout.write("%s\t%s\t%s\t%s\n" % (chrom2, line[2], line[2], "\t".join(line[3:])))
        else:
            sys.stdout.write(s)
