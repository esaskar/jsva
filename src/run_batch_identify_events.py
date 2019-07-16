#!/usr/bin/env python2

import sys, subprocess, os

import sbrun

LIST="jsva_cnv.file_matching"
ODIR="jsva.complex.end_matching.no_multisample.r100k"

try:
    os.makedirs(ODIR)
except:
    pass

for s in open(LIST):
    sample, jsvafn, cnvfn = s.strip().split("\t")
    print sample, jsvafn, cnvfn
    outfn = "%s/%s.out" % (ODIR, sample)
    errfn = "%s/%s.err" % (ODIR, sample)
    jobfn = "%s/%s.job" % (ODIR, sample)
    sbrun.run_cmd("src/identify_events.py -i %s -c %s -o %s/%s -t 7 -n -z -m 5 -r 100000 --end-matching --omit-multisample" % (jsvafn, cnvfn, ODIR, sample), vmem = 3000, out_fn = outfn, err_fn = errfn, job_fn = jobfn)
#    sys.exit(0)