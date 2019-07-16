#!/usr/bin/env python

import sys

f = sys.stdin
o = sys.stdout
for s in f:
    if s.startswith("#"):
        continue
    o.write(s)
    break

svids = {}
np = n = 0
for s in f:
    svid = s.strip().split()[4]
    if svid in svids:
        o.write(svids[svid])
        o.write(s)
        del svids[svid]
        np += 1
    else:
        svids[svid] = s
        n += 1

#print "%d pairs, %d total" % (np, n)
