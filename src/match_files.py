#!/usr/bin/env python2
"""
This script attempts to match files in two directories A and B by
first checking whether the basename (without a suffix) of file in A occurs as a substring
in filename in B (e.g., "c112.gz" -> "c112" occurs in "Fam_c112_1_5980_T_LP6005743-DNA_D03.somatic.cnv.bed.gz"),
and then vice versa.

Author: Esa Pitkanen (esa.pitkanen@helsinki.fi)
"""

import sys, os, argparse

def remove_suffix(fn):
    k = fn.rfind(".")
    if k != -1:
        return fn[:k]
    else:
        return fn

def __main(first, second, output, suffix):
    fns1 = filter(lambda x: suffix == None or x.endswith(suffix), os.listdir(first))
    fns2 = filter(lambda x: suffix == None or x.endswith(suffix), os.listdir(second))
    if output:
        o = open(output, "w")
    else:
        o = sys.stdout

    matches = set()
    no_match1 = set(fns1)
    no_match2 = set(fns2)
        
    n = 0
    for fn1 in fns1:
        s1 = remove_suffix(fn1)
        for fn2 in fns2:
            if os.path.basename(s1) in fn2:
                o.write("%s/%s\t%s/%s\n" % (first, fn1, second, fn2))
                n += 1
                matches.add(fn1)
                no_match1.discard(fn1)
                no_match2.discard(fn2)
                
    for fn2 in fns2:
        s2 = remove_suffix(fn2)
        for fn1 in fns1:
            if fn1 in matches:
                continue
            if os.path.basename(s2) in fn1:
                o.write("%s/%s\t%s/%s\n" % (first, fn1, second, fn2))
                n += 1
                no_match1.discard(fn1)
                no_match2.discard(fn2)

    sys.stderr.write("%d files in %s, %d files in %s, %d matches\n" % (len(fns1), first, len(fns2), second, n))
    for fn in no_match1:
        sys.stderr.write("No match\t%s\t%s\n" % (fn, first))
    for fn in no_match2:
        sys.stderr.write("No match\t%s\t%s\n" % (fn, second))

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-1", "--first", required = True)
    p.add_argument("-2", "--second", required = True)
    p.add_argument("-o", "--output")
    p.add_argument("-s", "--suffix", help = "Files must have this suffix to match")
    args = p.parse_args()
    __main(args.first, args.second, args.output, args.suffix)
    
