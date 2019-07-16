#!/bin/bash

DATA=complex.k10000.opposite_breakends.txt

# pie graph of event types
#tail -n +2 $DATA | cut -f 2 | sort | uniq -c | awk '{print $2,$1}' > tmp/complex.types

# stacked barplot annotations vs event types

Rscript src/plot_complex.R
