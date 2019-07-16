#!/usr/bin/env python

import sys, subprocess, os, argparse

p = argparse.ArgumentParser()
p.add_argument("-f", "--flank", default = 10000, type = int)
args = p.parse_args()

SAMPLES = "data/samples_with_cohesin_mutation"
#COHESIN_SITES = "data/CTCF_RAD21_intersect.bed"
#COHESIN_SITES = "data/occupied_cohesin_sites_with_mutation.bed"
COHESIN_SITES = "data/CTCF_RAD21_intersect_mutated.bed"
SITE_MUTATIONS = "data/somatic_mutations_occupied_cohesin.bed"
BREAKPOINTS = "out/somatic_breakpoints/sorted/"

o1 = open("breakpoints_at_mutated_sites.txt", "w")
o2 = open("breakpoints_at_nonmutated_sites.txt", "w")
obpcm = open("breakpoint_coords_at_mutated_sites.txt", "w")
obpcn = open("breakpoint_coords_at_nonmutated_sites.txt", "w")

# 0. read samples
samples = map(lambda x: x.strip(), open(SAMPLES).readlines())
samples.sort()
print "%d samples from %s" % (len(samples), SAMPLES)

# map samples to somatic breakpoint files
sample_breakpoints = {}
for s in samples:
    fn = "%s/%s.gz" % (BREAKPOINTS, s)
    if os.path.exists(fn):
        sample_breakpoints[s] = fn
        continue
    fn = "%s/%s_1.gz" % (BREAKPOINTS, s)
    if os.path.exists(fn):
        sample_breakpoints[s] = fn
        continue
    print "Can't find", s
    exit()

class Site:
    def add_sample(self, sample):
        self.samples.add(sample)
    def __init__(self, chrom, start, stop, orient):
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.orient = orient
        self.samples = set()
    def __str__(self):
        return "%s:%d-%d:%s" % (self.chrom, self.start, self.stop, self.orient)

# 1. Collect occupied cohesin sites
pos_to_site = {}
for s in open(COHESIN_SITES):
    v = s.strip().split("\t")
    chrom, start, stop, tag, value, orientation = v
    start, stop = int(start), int(stop)
    site = Site(chrom, start, stop, orientation)
    key = "%s:%d-%d" % (chrom, start, stop)
    pos_to_site[key] = site

print "%d occupied cohesin sites from %s" % (len(pos_to_site), COHESIN_SITES)

# 2. Identify samples+sites with mutation
n = 0
for s in open(SITE_MUTATIONS):
    if s.startswith("#"):
        continue
    chrom, start, end, sample = s.strip().split("\t")
    start, end = int(start), int(end)
    key = "%s:%d-%d" % (chrom, start, end)
    if key not in pos_to_site:
        print "site", key, "not found"
        exit()
    pos_to_site[key].add_sample(sample)
    n += 1
print "%d mutations in sites from %s" % (n, SITE_MUTATIONS)
#print "%d/%d mutated sites" % ()

no_mut = mut = 0
for k in pos_to_site:
    if len(pos_to_site[k].samples) == 0:
        no_mut += 1
    else:
        mut += 1

print "%d/%d sites mutated in >0 samples" % (mut, no_mut + mut)

# 3. Count SV/CNV breakpoints near samples+sites

def get_breakpoints_in_neighborhood(site, sample, flank):
    chrom = site.chrom
    start = site.start - flank
    end = site.stop + flank
    s = subprocess.check_output("tabix %s %s:%d-%d" % (sample_breakpoints[sample], chrom, start, end), shell = True)
    v = s.strip().split("\n")
    n = len(s.strip().split("\n")) - 1
    assert(n >= 0)
    if len(v) <= 1:
        s = None
    else:
        s = []
        for vals in v:
            chrom, pos, svtype, reads, qual = vals.split("\t")
            pos = int(pos)
            # tpos > 0 iff breakpoint to 3' of site start
            if site.orient == "+":
                tpos = pos - site.start
            elif site.orient == "-":
                tpos = site.stop - pos
            s.append("\t".join(map(str, [chrom, pos, tpos, svtype, reads, qual, sample, site.orient])))
        s = "%s\n" % ("\n".join(s))
    return n, s

mutated_site_breakpoints = []
nonmutated_site_breakpoints = []

keys = pos_to_site.keys()
keys.sort()
for i, k in enumerate(keys):
    site = pos_to_site[k]
    nm = mm = 0
    totaln = totalm = 0
    for sample in samples:
        #print sample, sample in site.samples
        n, bps = get_breakpoints_in_neighborhood(site, sample, args.flank)
 #       print sample, sample in site.samples, "%d breakpoints" % (n)
        if sample in site.samples:
            mutated_site_breakpoints.append(n)
            if bps != None:
                obpcm.write(bps)
            mm += n
            totalm += 1
        else:
            nonmutated_site_breakpoints.append(n)
            if bps != None:
                obpcn.write(bps)
            nm += n
            totaln += 1
    print "%d/%d" % (i, len(keys)), k, nm, totaln, mm, totalm
#    if i > 1000:
#        break

o1.write("\n".join(map(str, mutated_site_breakpoints)))
o2.write("\n".join(map(str, nonmutated_site_breakpoints)))
