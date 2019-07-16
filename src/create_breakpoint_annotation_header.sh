#!/bin/bash

BEDS="UCSC_Genes_knownGene.hg19.nochr.bed UCSC_Pfam_ucscGenePfam.nameOnly.hg19.bed ENCODE.bed RepeatMasker.hg19.bed DGV.hg19.bed Lovo_sorted.bed dukeExcludeRegions_hg19.nochr.bed HiSeqDepth.hg19.0.05.bed Mappability.hg19.bed wgEncodeCrgMapabilityAlign100mer.threeormore.merged.flank_r99_hg19.nochr.bed"

FN=breakpoint_annotation_header_temp.txt
echo "#Header" > $FN
for b in $BEDS; do 
    echo -en "$b\t" >>$FN
done
echo >>$FN

for b in $BEDS; do 
    echo -en "$(head -n 1 BED/$b)\t" >>$FN
done

echo Edit $FN
