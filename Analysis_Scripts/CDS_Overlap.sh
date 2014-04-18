#!/bin/bash

#   A script to generate the CDS file against which we annotate SNPs
#   as coding or noncoding. Takes two arguments: 
#       1) The reduced GFF file that contains only CDS annotations
#       2) The BED file that contains all the SNPs 
#   This script will produce a modified GFF file that is suitable for 
#   the next script in the pipeline (Filter_CDS.py)

#   Check if we have the correct number of arguments
if [ $# != 2 ]
    then
    echo "Usage:"
    echo "`basename $0` GFF BED > Mod_GFF"
    echo "Will produce a modified GFF that contains only the CDS annotations"
    echo "that overlap variants listed in the BED file. This will then be used"
    echo "by the next script in the pipeline to annotate SNPs with amino acid"
    echo "and codon states."
    echo ""
    echo "Requires bedtools."
    exit 1
fi

#   Check for bedtools
command -v bedtools > /dev/null 2>&1 || { echo >&2 "This script requires bedtools."; exit 1; }

#   Set the arguments to have meaningful names
GFF=$1
BED=$2

bedtools intersect -a $BED -b $GFF -wo | awk '{print $1 "\t" $9 "\t" $10 "\t" $5 "\t" $12 "\t" $13 }'
