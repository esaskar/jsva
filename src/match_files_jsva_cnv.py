#!/usr/bin/env python2
"""
This script attempts to match files in two directories A and B by
first checking whether the basename (without a suffix) of file in A occurs as a substring
in filename in B (e.g., "c112.gz" -> "c112" occurs in "Fam_c112_1_5980_T_LP6005743-DNA_D03.somatic.cnv.bed.gz"),
and then vice versa.

This version of match_files.py is for matching VarScan CRC files against JSVA files.

Author: Esa Pitkanen (esa.pitkanen@helsinki.fi)
"""

import sys, os, argparse

def remove_suffix(fn):
    if fn.startswith("Fam_"):
        fn = fn[4:]
    k = fn.rfind(".")
    if k != -1:
        fn = fn[:k]
    fn = fn.split("_")[0]
    return fn

def __main(first, second, output, suffix):
    fns1 = filter(lambda x: suffix == None or x.endswith(suffix), os.listdir(first))
    fns2 = filter(lambda x: suffix == None or x.endswith(suffix), os.listdir(second))
    if output:
        o = open(output, "w")
    else:
        o = sys.stdout

    matches = {}
    no_match1 = set(fns1)
    no_match2 = set(fns2)
        
    n = 0
    for fn1 in fns1:
        s1 = remove_suffix(os.path.basename(fn1))
        for fn2 in fns2:
            s2 = remove_suffix(os.path.basename(fn2))
            if s1 == s2:
                o.write("%s\t%s/%s\t%s/%s\n" % (s1, first, fn1, second, fn2))
                n += 1
                if fn1 not in matches:
                    matches[fn1] = [fn2]
                else:
                    matches[fn1].append(fn2)
                if fn2 not in matches:
                    matches[fn2] = [fn1]
                else:
                    matches[fn2].append(fn1)
                no_match1.discard(fn1)
                no_match2.discard(fn2)
                
    for fn2 in fns2:
        s2 = remove_suffix(os.path.basename(fn2))
        for fn1 in fns1:
            if fn1 in matches:
                continue
            s1 = remove_suffix(os.path.basename(fn1))
            if s1 == s2:
                o.write("%s/%s\t%s/%s\n" % (s1, first, fn1, second, fn2))
                n += 1
                no_match1.discard(fn1)
                no_match2.discard(fn2)
                if fn1 not in matches:
                    matches[fn1] = [fn2]
                else:
                    matches[fn1].append(fn2)
                if fn2 not in matches:
                    matches[fn2] = [fn1]
                else:
                    matches[fn2].append(fn1)

    for fn in no_match1:
        sys.stderr.write("No match\t%s\t%s\n" % (fn, first))
    for fn in no_match2:
        sys.stderr.write("No match\t%s\t%s\n" % (fn, second))

    k = matches.keys()
    k.sort()
    for fn in k:
        if len(matches[fn]) > 1:
            sys.stderr.write("Multiple\t%s\t%s\n" % (fn, ",".join(matches[fn])))

    sys.stderr.write("%d files in %s, %d files in %s, %d matches\n" % (len(fns1), first, len(fns2), second, n))
        
if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-1", "--first", required = True)
    p.add_argument("-2", "--second", required = True)
    p.add_argument("-o", "--output")
    p.add_argument("-s", "--suffix", help = "Files must have this suffix to match")
    args = p.parse_args()
    __main(args.first, args.second, args.output, args.suffix)
    
