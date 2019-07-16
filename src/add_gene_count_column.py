#!/usr/bin/env python

import sys

f = sys.stdin
o = sys.stdout

prev_gene = None
data = []
samples = set()
for s in f:
    v = s.strip().split("\t")
    gene = v[23]
    if prev_gene == None:
        prev_gene = gene
    if gene == prev_gene:
        data.append(v)
        samples.add(v[5])
    else:
        n = len(samples)
        for d in data:
            o.write("%s\t%s\n" % ("\t".join(d), n))
        data = []
        samples = set()
        prev_gene = gene

for d in data:
    o.write("%s\t%s\n" % ("\t".join(d), n))
