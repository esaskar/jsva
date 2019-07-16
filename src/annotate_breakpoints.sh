#!/bin/bash

set -euo pipefail

BEDDIR="$1"
BEDS="sorted_euL1db_full_L1.bed" 

# jsva library
SORTED="simple_events_by_breakpoint_sorted.tsv.gz"

first=$(tempfile)
echo "Copying $SORTED -> $first"
cp $SORTED $first

for bed in $BEDS; do
    bfn=$BEDDIR/$bed
    second=$(tempfile)
    echo "$bfn: $first -> $second"
    bedtools closest -a $first -b $bfn -D ref -t first > $second
    if [[ $? != 0 ]]; then
	exit $?;
    fi
    rm $first
    first=$second
done

cp $second $RESULT
