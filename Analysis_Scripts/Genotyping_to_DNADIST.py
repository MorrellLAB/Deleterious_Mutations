#!/usr/bin/env python

#   A script to go from the tabular genotyping format to the format expected
#   for input to dnadist in PHYLIP. Will insert N for heterozygous SNPs and
#   Ns for missing data

import sys

seq = ''
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        tmp = line.strip().split('\t')
        #   The first line contains the sample
        if index == 0:
            sample = tmp[1]
        else:
            #   Ugly, but the only good way to do it
            call = list(set(tmp[1]))
            #   If the length of call is greater than 1, or there is a -
            #   then set it to N
            if len(call) > 1 or '-' in call:
                seq = seq + 'N'
            else:
                seq = seq + call[0]

#   Then print the sequence out
print sample + '\t' + seq
