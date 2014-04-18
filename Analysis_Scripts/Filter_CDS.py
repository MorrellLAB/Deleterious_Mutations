#!/usr/bin/env python

#   A very simple intermediary script to remove duplicate CDS entries
#   from the CDS_Overlap.sh script. Bedtools sometimes reports multiple
#   hits from the calculation

import sys

snps = []
contents = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        tmp = line.split('\t')
        snp_id = tmp[3]
        if snp_id in snps:
            continue
        else:
            contents.append(line)
        snps.append(snp_id)
print ''.join(contents)
