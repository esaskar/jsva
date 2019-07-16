#!/usr/bin/env python

import sys, os, argparse, gzip, errno

p = argparse.ArgumentParser()
p.add_argument("-i", "--input")
p.add_argument("-o", "--output")
p.add_argument("-e", "--header")
p.add_argument("-s", "--show-header", action = "store_true")
p.add_argument("-z", "--zero", action = "append", default = [], help = "Reject if x != 0")
p.add_argument("-n", "--nonzero", action = "append", default = [], help = "Reject if x == 0")
p.add_argument("-g", "--greater", action = "append", default = [], help = "Reject if x <= v")
p.add_argument("-l", "--lesser", action = "append", default = [], help = "Reject if x >= v")
p.add_argument("-m", "--match", action = "append", default = [], help = "Reject if x != v")
p.add_argument("-x", "--mismatch", action = "append", default = [], help = "Reject if x == v")
p.add_argument("-c", "--proximal", action = "append", default = [], help = "Reject if abs(x) > v")
p.add_argument("-f", "--distal", action = "append", default = [], help = "Reject if abs(x) < v")
args = p.parse_args()

fn = args.input
if fn == None or fn == "-":
    f = sys.stdin
elif fn.endswith(".gz") == False:
    f = open(fn)
else:
    f = gzip.open(fn)

#if args.header:
#    headerfn = args.header
#else:
#    headerfn = fn.rstrip(".tsv.gz") + ".header.txt"    
#if os.path.exists(headerfn) == False:
#    sys.stderr.write("Header file %s missing\n" % (headerfn))
#    sys.exit(1)

args.greater = map(lambda x: (x[0], float(x[1])), map(lambda x: x.split(","), args.greater))
args.lesser = map(lambda x: (x[0], float(x[1])), map(lambda x: x.split(","), args.lesser))
args.proximal = map(lambda x: (x[0], float(x[1])), map(lambda x: x.split(","), args.proximal))
args.distal = map(lambda x: (x[0], float(x[1])), map(lambda x: x.split(","), args.distal))
args.match = map(lambda x: x.split(","), args.match)
args.mismatch = map(lambda x: x.split(","), args.mismatch)

#header = open(headerfn).read().split()
s = f.readline()
while s.startswith("#"):
    s = f.readline()
header = s.strip().split()

if args.show_header:
    print "\n".join(header)
    sys.exit(0)

if args.output:
    o = open(args.output, "w")
else:
    o = sys.stdout

columns = {}
cols = []
cols.extend(args.zero)
cols.extend(args.nonzero)
cols.extend(map(lambda x: x[0], args.greater))
cols.extend(map(lambda x: x[0], args.lesser))
cols.extend(map(lambda x: x[0], args.match))
cols.extend(map(lambda x: x[0], args.mismatch))
cols.extend(map(lambda x: x[0], args.distal))
cols.extend(map(lambda x: x[0], args.proximal))
sys.stderr.write("%d columns specified\n" % (len(cols)))
if len(cols) == 0:
    sys.stderr.write("Specify at least one column\n")
    sys.exit(2)
for i, k in enumerate(header):
    columns[k] = i
for c in cols:
    if c not in columns:
        sys.stderr.write("Column %s not found in header\n" % (c))
        sys.exit(1)

o.write("%s\n" % ("\t".join(header)))

n = num_passed = 0
try:
    for s in f:
        v = s.strip().split("\t")
        n += 1
        passed = True
        for c, x in args.match:
            if v[columns[c]] != x:
                passed = False
        for c, x in args.mismatch:
            if v[columns[c]] == x:
                passed = False
        for c in args.zero:
            if int(v[columns[c]]) != 0:
                passed = False
        for c in args.nonzero:
            if int(v[columns[c]]) == 0:
                passed = False
        for c, x in args.greater:
            if float(v[columns[c]]) <= x:
                passed = False
        for c, x in args.lesser:
            if float(v[columns[c]]) >= x:
                passed = False
        for c, x in args.proximal:
            if abs(float(v[columns[c]])) > x:
                passed = False
        for c, x in args.distal:
            if abs(float(v[columns[c]])) < x:
                passed = False
        if passed:
            o.write(s)
            num_passed += 1
except IOError, e:
    if e.errno == errno.EPIPE:
        pass
    else:
        raise

sys.stderr.write("%d/%d (%.1f%%) passed breakpoints\n" % (num_passed, n, 100.0 * num_passed / n))
