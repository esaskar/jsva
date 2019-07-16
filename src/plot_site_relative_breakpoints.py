#!/usr/bin/env python

import sys, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

f1 = open(sys.argv[1])  # breakpoint coords for mutated samples
f2 = open(sys.argv[2])  # breakpoint coords for nonmutated samples

minv = sys.maxint
maxv = -sys.maxint

def parse_bp(f):
    global minv, maxv
    bps = []
    svbp = {"del" : [], "dup" : [], "inv" : [], "tr" : []}
    for s in f:
        v = s.strip().split("\t")
        relpos = int(v[2])
        svtype = v[3]
        reads, qual = int(v[4]), int(v[5])
        minv = min(minv, relpos)        
        maxv = max(maxv, relpos)
        bps.append((relpos, svtype))
        svbp[svtype].append(relpos)
    return bps, svbp

muts, svmuts = parse_bp(f1)
nonmuts, svnonmuts = parse_bp(f2)

def get_color(svtype):
    if svtype == "del":
        c = "red"
    elif svtype == "dup":
        c = "green"
    elif svtype == "inv":
        c = "blue"
    elif svtype == "tr":
        c = "orange"
    else:
        c = "black"
    return c

from scipy.stats import gaussian_kde

def plot_density(ax, x, col, divisor):
    if len(x) == 0:
        return
    xs = np.linspace(minv, maxv, 1000)
    d = gaussian_kde(x, bw_method = .01)
#    d.covariance_factor = lambda : .1
#    d._compute_covariance()
    ys = d(xs) * len(x) / divisor
#    ys = d(xs)
    ax.plot(xs, ys, color = col, lw = 0)
    ax.fill_between(xs, ys, alpha = .5, zorder = 5, antialiased = True, color = col)
    ax.plot(x, np.zeros(len(x)), ms = 5, color = "black", marker = "|")

def plot_sub(ax, data, div):
    plot_density(ax, data["dup"], "green", div)
    plot_density(ax, data["del"], "red", div)
    plot_density(ax, data["inv"], "blue", div)
    plot_density(ax, data["tr"], "orange", div)

fig = plt.figure()
ax1 = fig.add_subplot(2, 1, 1)
plot_sub(ax1, svmuts, 1.0)
ax1.set_xlim(minv, maxv)

ax2 = fig.add_subplot(2, 1, 2)
plot_sub(ax2, svnonmuts, 1.0)
ax2.set_xlim(minv, maxv)

plt.savefig("site_breakpoints.png", dpi = 300, bbox_inches = "tight")

fig = plt.figure()
plot_density(plt, svmuts["dup"] + svmuts["del"] + svmuts["inv"] + svmuts["tr"], "red", 984)
plot_density(plt, svnonmuts["dup"] + svnonmuts["del"] + svnonmuts["inv"] + svnonmuts["tr"], "blue", 77324)
plt.gca().set_xlim(minv, maxv)
plt.savefig("site_breakpoints_combined.png", dpi = 300, bbox_inches = "tight")

fig = plt.figure()
bins = np.linspace(minv, maxv, 15)
nm = np.array(svnonmuts["dup"] + svnonmuts["del"] + svnonmuts["inv"] + svnonmuts["tr"]) #/ 77324.0
mm = np.array(svmuts["dup"] + svmuts["del"] + svmuts["inv"] + svmuts["tr"]) #/ 984.0
n2, bins2, patches2 = plt.hist(nm, bins = bins, facecolor = "blue", lw = 0, alpha = 0.5, normed = 1)
n1, bins1, patches1 = plt.hist(mm, bins = bins, facecolor = "red", lw = 0, alpha = 0.5, normed = 1)
plt.gca().set_xlim(minv, maxv)
plt.savefig("site_breakpoints_combined_histogram.png", dpi = 300, bbox_inches = "tight")
