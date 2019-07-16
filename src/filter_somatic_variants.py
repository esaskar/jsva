#!/usr/bin/env python

import sys, argparse, gzip, re

NUM_FIRST_COLUMNS = 7  # columns before tumor types
NUM_LAST_COLUMNS = 11  # columns after individuals

CALLFN = "WGS_SV_simple_calls_20140513.csv.gz"

p = argparse.ArgumentParser()
p.add_argument("-i", "--input", help = "SV call file [%(default)s]", default = CALLFN)
p.add_argument("-o", "--output")
p.add_argument("-t", "--tumortype")
p.add_argument("-n", "--normal-freq", help = "Max population frequency to allow in germline [%(default)s]", default = 0.01, type = float)
args = p.parse_args()

if args.input.endswith(".gz") or args.input.endswith(".gzip"):
    f = gzip.open(args.input)
else:
    f = open(args.input)

if args.output != None:
    o = open(args.output, "w")
else:
    o = sys.stdout

if args.tumortype != None:
    args.tumortype = set(args.tumortype.split(","))

rett = re.compile("(\w+)/(\d+)")

hdr = f.readline().strip().split("\t")
ng_total = int(hdr[5].split("/")[1])
ns_total = int(hdr[6].split("/")[1])
tumortypes = {}
selected_tumor_types = set()
individuals = {}
for i, v in enumerate(hdr):
    if i >= NUM_FIRST_COLUMNS:
        tt = rett.findall(v)
        if len(tt) > 0:
            # tumor types
            tt = tt[0]
            tname, ct = tt[0].split("_")
            tumortypes[i] = int(tt[1])
            if args.tumortype != None and tname in args.tumortype and ct == "S":
                selected_tumor_types.add(i)
                print tname, i
        elif i < len(hdr) - NUM_LAST_COLUMNS:
            # individuals
            assert(v not in individuals)
            individuals[v] = i

print "\t".join(hdr)

sys.stderr.write("Max germline freq %f -> max %d germline hits\n" % (args.normal_freq, args.normal_freq * ng_total))

for s in f:
    v = s.strip().split("\t")
    chrom, bp1, bp2, length, svtype, ngermline, nsomatic = v[:7]
    ngermline, nsomatic = int(ngermline), int(nsomatic)
    gfreq = 1.0 * ngermline / ng_total
    if gfreq < args.normal_freq:
        ok = False
        if selected_tumor_types != None:
            for i in selected_tumor_types:
                nt = int(v[i])
                if nt > 0:
                    print v[i]
                    ok = True
                    break
        else:
            ok = True
        if ok:
            print s.strip()
