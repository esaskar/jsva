#!/usr/bin/env python2
"""
Annotated structural variant breakpoint
"""

import gzip

class Breakpoint:
    SV_DEL = 0
    SV_INV = 1
    SV_DUP = 2
    SV_TR = 3
    SV_TYPE = {"del" : SV_DEL,
               "inv" : SV_INV,
               "dup" : SV_DUP,
               "tr" : SV_TR}
    SV_TYPE_NAME = {SV_DEL : "del",
                    SV_INV : "inv",
                    SV_DUP : "dup",
                    SV_TR : "tr"}
    BREAKEND_5P = "-"  #       |---
    BREAKEND_3P = "+"  #  ---|

    def __init__(self, s):
        v = s.split("\t")
        self.chrom, self.pos = v[0], int(v[1])
        self.length = v[3]
        if self.length != "-":
            self.length = int(self.length)
        self.svid = v[4]
        self.svtype = Breakpoint.SV_TYPE[v[5]]
        self.ngermline = int(v[6])
        self.nsomatic = int(v[7])
        self.data = v
        self.orient = None  # resolved while parsing data in BreakpointSet
        self.samples = None # resolved  "

    def tumor_read_counts(self):
        tumor_reads = []
        for sample in self.calls:
            call_type, normal, tumor = self.calls[sample]
            if tumor != "-":
                tumor_reads.append(int(tumor.split(":")[0]))
        return tumor_reads

    def __repr__(self):
        return "%s:%s:%s:%s:%s" % (self.chrom, self.pos, self.svid, Breakpoint.SV_TYPE_NAME[self.svtype], self.orient)
    def __hash__(self):
        return hash("%s:%s" % (self.svid, self.pos))

def nearest_breakpoint(bp, bp1chr, bp1pos, bp2chr, bp2pos):
    # Return 1 if bp1, 2 if bp2 is closest to bp
    if bp.chrom != bp1chr:
        if bp.chrom != bp2chr:
            return None
        else:
            return 2
    else:
        if bp.chrom != bp2chr:
            return 1
        else:
            if abs(bp.pos - bp1pos) < (bp.pos - bp2pos):
                return 1
            else:
                return 2

class BreakpointSet:
    IX_BP1_CHR = "BP1Chr"
    IX_BP1_POS = "BP1Mean"
    IX_BP1_OR = "BP1Or"
    IX_BP2_CHR = "BP2Chr"
    IX_BP2_POS = "BP2Mean"
    IX_BP2_OR = "BP2Or"
    IX_LENGTH_MEAN = "LengthMean"

    FIRST_SAMPLE_SET = 8

    TAG_GERMLINE = "G"
    TAG_SOMATIC = "S"

    def __init__(self, df, header_fn = None, 
                 max_germline_af = 1.0, 
                 min_normal_reads = 0, min_tumor_reads = 0):
        self.max_germline_af = max_germline_af
        self.min_normal_reads = min_normal_reads
        self.min_tumor_reads = min_tumor_reads
        self.bpf = df
        if header_fn:
            hdrf = open(header_fn)
        else:
            hdrf = self.bpf
        self.header = self.__parse_header(hdrf)
        self.header_ix = {}
        for i, h in enumerate(self.header):
            self.header_ix[h] = i
        self.num_germline = int(self.header[6].split("/")[1])
        self.num_somatic = int(self.header[7].split("/")[1])

        # iterate over sample sets in header and parse number of samples in each set
        self.sample_sets = []
        self.n_samples = {}
        self.sample_of_set = {}
        self.samples = []
        for i in range(BreakpointSet.FIRST_SAMPLE_SET, len(self.header)):
            if "/" not in self.header[i]:
                # first sample encountered
                self.first_sample = i
                break
            sample_set, n = self.header[i].split("/")
            set_name, germline_somatic = sample_set.split("_")
            n = int(n)
            if (set_name, n) not in self.sample_sets:
                self.sample_sets.append((set_name, n))
                self.sample_of_set[set_name] = []

            assert(germline_somatic == BreakpointSet.TAG_GERMLINE or germline_somatic == BreakpointSet.TAG_SOMATIC)
            if set_name not in self.n_samples:
                self.n_samples[set_name] = {BreakpointSet.TAG_GERMLINE : 0, BreakpointSet.TAG_SOMATIC : 0}
            self.n_samples[set_name][germline_somatic] = n

        # iterate over samples
        # TODO: current data version does not enumerate samples in the same order as sample sets :(
        self.last_sample = self.header_ix[BreakpointSet.IX_LENGTH_MEAN] - 1
        i = 0
        c = self.first_sample
        for set_name, n in self.sample_sets:
            while i < n:
                self.samples.append(self.header[c])
                self.sample_of_set[set_name].append(self.header[c])
                c += 1
                i += 1
            i = 0

    def get_column(self, name):
        return self.header_ix[name]

    def __parse_header(self, f):
        while 1:
            s = f.readline()
            if s.startswith("#") == False:
                return s.strip().split()
        raise Exception("EOF while scanning for header")

    def __iter__(self):
        return self

    def next(self):
        while 1:
            s = self.bpf.readline().strip()
            if s == "":
                raise StopIteration
            bp = Breakpoint(s)
            bp.bp1chr = bp.data[self.header_ix[BreakpointSet.IX_BP1_CHR]]
            bp.bp1pos = float(bp.data[self.header_ix[BreakpointSet.IX_BP1_POS]])
            bp.bp1or = bp.data[self.header_ix[BreakpointSet.IX_BP1_OR]]
            bp.bp2chr = bp.data[self.header_ix[BreakpointSet.IX_BP2_CHR]]
            bp.bp2pos = float(bp.data[self.header_ix[BreakpointSet.IX_BP2_POS]])
            bp.bp2or = bp.data[self.header_ix[BreakpointSet.IX_BP2_OR]]

            nearest_bp = nearest_breakpoint(bp, bp.bp1chr, bp.bp1pos, bp.bp2chr, bp.bp2pos)
            if nearest_bp == None:
                sys.stderr.write("Warning: can't resolve orientation for breakpoint %s\n" % (bp))
            else:
                if nearest_bp == 1:
                    bp.orient = bp.bp1or
                else:
                    bp.orient = bp.bp2or

            bp.germline_samples = set()
            bp.somatic_samples = set()
            bp.calls = {}
            n_germline_calls = 0
            n_calls = 0
            for i in range(self.first_sample, self.last_sample):
                if bp.data[i] != "-":
                    n_calls += 1
                    si = i - self.first_sample
                    call_type, normal, tumor = bp.data[i].split("/")
                    if normal != "-":
                        n_germline_calls += 1
                        n_reads, n_qual, n_range = normal.split(":")
                        n_reads, n_qual = int(n_reads), int(n_qual)
                    else:
                        n_reads = n_qual = 0
                    if tumor != "-":
                        t_reads, t_qual, t_range = tumor.split(":")
                        t_reads, t_qual = int(t_reads), int(t_qual)
                    else:
                        t_reads = t_qual = 0
                    if n_reads < self.min_normal_reads or t_reads < self.min_tumor_reads:
                        continue
                    sample = self.samples[si]
                    if normal != "-":
                        bp.germline_samples.add(sample)
                    elif tumor != "-":
                        bp.somatic_samples.add(sample)
                    bp.calls[sample] = (call_type, normal, tumor)
    #                print self.samples[si], bp.data[i]

            if n_calls == 0 or 1.0 * n_germline_calls / n_calls > self.max_germline_af:
                continue
            return bp

if __name__ == "__main__":
    import sys
    bps = BreakpointSet(open(sys.argv[1]))
    for s in bps:
        print s
