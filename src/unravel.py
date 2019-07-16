#!/usr/bin/env python

import sys, argparse, os, gzip, errno, datetime

class SampleBreakpoint:
    def __init__(self, sstr = None):
        v = sstr.strip().split("\t")
        self.chrom, self.start, self.end, self.length, self.svid, self.svtype, self.or1, self.or2, self.sample, self.n_germline, self.n_somatic, self.nreads, self.nqual, self.nbp1, self.nbp2, self.treads, self.tqual, self.tbp1, self.tbp2, self.mbp1chrom, self.mbp1pos, self.mbp2chrom, self.mbp2pos = v[:23]
        self.n_somatic = int(self.n_somatic)
        self.n_germline = int(self.n_germline)
        self.nreads = int(self.nreads)
        self.treads = int(self.treads)
        if len(v) > 27:
            self.gene = (v[26], int(v[27]))
        else:
            self.gene = (None, None)
        if len(v) > 35:
            self.encode = (v[34], int(v[35]))
        else:
            self.encode = (None, None)
        if len(v) > 44:
            self.repeat = (v[43], int(v[44]))
        else:
            self.repeat = (None, None)

    def get_label(self, max_dist):
        if self.gene[1] <= max_dist:
            return self.gene[0]
        elif self.encode[1] <= max_dist:
            return self.__conv_encode_label(self.encode[0])
        else:
            return None

    def __conv_encode_label(self, x):
        return "_".join(x.split("_")[1:2]) # 7_Weak_Enhancer -> Weak_Enhancer

    def __repr__(self):
        return "%s:%s:%s" % (self.svtype, self.chrom, self.start)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input")
    p.add_argument("-o", "--output")
    p.add_argument("-e", "--header")
    p.add_argument("-r", "--reads", help = "Minimum number of reads in tumor", default = 0, type = int)
    p.add_argument("-q", "--qual", help = "Minimum quality in tumor", default = 0, type = int)
    p.add_argument("--normal-reads", help = "Minimum number of reads in normal", default = 0, type = int)
    p.add_argument("--normal-qual", help = "Minimum quality in tumor", default = 0, type = int)
    args = p.parse_args()

    # if args.header:
    #     headerfn = args.header
    # else:
    #     headerfn = args.input[:-len(".tsv.gz")] + ".header.txt"    
    # if os.path.exists(headerfn) == False:
    #     sys.stderr.write("Header file %s missing\n" % (headerfn))
    #     sys.exit(1)
    #print header[first_sample_ix], header[last_sample_ix]

    if args.input and args.input != "-":
        if args.input.endswith(".gz"):
            f = gzip.open(args.input)
        else:
            f = open(args.input)
    else:
        f = sys.stdin

    if args.output:
        o = open(args.output, "w")
    else:
        o = sys.stdout

#    header = open(headerfn).read().split()
    while 1:
        s = f.readline()
        if s.startswith("#"):
            continue
        else:
            header = s.strip().split()
            break
    columns = {}
    first_sample_ix = last_sample_ix = None
    for i, k in enumerate(header):
        columns[k] = i
        if "/" not in k and i > 9 and first_sample_ix == None:
            first_sample_ix = i
        if k == "LengthMean":
            last_sample_ix = i - 1
    assert(first_sample_ix != None)
    assert(last_sample_ix != None)

    n_germline = int(header[10].split("/")[1])
    n_somatic = int(header[11].split("/")[1])

    def parse_call(call):
        if call == "-":
            reads = qual = 0
            bp1 = bp2 = "-"
        else:
            reads, qual, bps = call.split(":")
            bp1, bp2 = bps.split("-")
        return int(reads), int(qual), bp1, bp2

    OR1IX = last_sample_ix + 2
    OR2IX = last_sample_ix + 7
    BP1CHR = last_sample_ix + 2
    BP1POS = last_sample_ix + 3
    BP2CHR = last_sample_ix + 6 
    BP2POS = last_sample_ix + 7
    FIRST_ANNOTATION = last_sample_ix + 10
    #DISTANCE = last_sample_ix + 14 # EP: this was not used in the script, removed

    while 1:
        s = f.readline()
        if s.startswith("#"):
            continue
        else:
            break

    o.write("#Output of \"%s\" on %s\n" % (" ".join(sys.argv), datetime.datetime.now()))
    try:
        o.write("Chrom\tStart\tEnd\tLength\tSVID\tSVType\tOr1\tOr2\tSample\tGermlineCalls/%d\tSomaticCalls/%d\tNReads\tNQual\tNBP1\tNBP2\tTReads\tTQual\tTBP1\tTBP2\tMergedBP1Chrom\tMergedBP1PosMean\tMergedBP2Chrom\tMergedBP2PosMean\n") # t%s\n") % (n_germline, n_somatic, "\t".join(header[FIRST_ANNOTATION:])))
        for s in f:
            if s.startswith("#"):
                continue
            v = s.strip().split("\t")
            or1, or2 = v[OR1IX], v[OR2IX]
            bp1chr, bp1pos, bp2chr, bp2pos = v[BP1CHR], v[BP1POS], v[BP2CHR], v[BP2POS]
#            distance = int(v[DISTANCE]) # EP: this was not used in the script, removed
            germline_calls = int(v[10])
            somatic_calls = int(v[11])
            type = v[5]
            for i in range(first_sample_ix, last_sample_ix):
                call = v[i]
                if call == "-":
                    continue
                ctype, normal, tumor= call.split("/")
                nreads, nqual, nbp1, nbp2 = parse_call(normal)
                treads, tqual, tbp1, tbp2 = parse_call(tumor)
                if treads < args.reads or tqual < args.qual or nreads < args.normal_reads or nqual < args.normal_qual:
                    continue
                #o.write("%s\t%s\t%s\n" % ("\t".join(v[0:8]), "\t".join(map(str, [header[i], germline_calls, somatic_calls, nreads, nqual, nbp1, nbp2, treads, tqual, tbp1, tbp2, bp1chr, bp1pos, bp2chr, bp2pos])))) #, "\t".join(v[FIRST_ANNOTATION:])))
                o.write("%s\t%s\n" % ("\t".join(v[0:8]), "\t".join(map(str, [header[i], germline_calls, somatic_calls, nreads, nqual, nbp1, nbp2, treads, tqual, tbp1, tbp2, bp1chr, bp1pos, bp2chr, bp2pos])))) #, "\t".join(v[FIRST_ANNOTATION:])))
    except IOError, e:
        if e.errno != errno.EPIPE:
            raise


