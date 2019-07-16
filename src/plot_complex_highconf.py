#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = 15, 5

df = pd.read_csv("complex.k10000.opposite_breakends.nodupsamples.highconf.txt", sep = " ", header=None, names=["UniqueSamples", "EventType"])

x = df.groupby(["EventType", "UniqueSamples"])["UniqueSamples"].count().unstack("EventType").fillna(0)
x["counts"] = x.sum(axis = 1)
x.sort(columns = "counts", inplace = True)
print "Median structural event count:", x["counts"].median()
x.drop("counts", inplace = True, axis = 1)

pal = ["#fd8d3c", "#31a354", "#de2d26", "#2b8cbe", "#fecc5c"]

fig = plt.figure(figsize = (12, 6), dpi = 300)
ax = fig.add_subplot(111)
x.plot(kind = "bar", stacked = "True", lw = 0, grid = False, color = pal)
plt.xlabel("Samples")
plt.ylabel("Structural events")
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off
fig.set_size_inches(18.5,10.5)
plt.savefig("highconf_samples.png", dpi = 300, bbox_inches = "tight")
