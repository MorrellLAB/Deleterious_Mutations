#!/usr/bin/env python
"""Calculate per-sample observered heterozygosity as the proportion of all
polymorphic sites that are heterozygous in a VCF."""

import sys

#   Empty dictionary to store the numbers of homozygous and heterozygous
#   sites
het = {}

nsites = 0
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            tmp = line.strip().split()
            samples = tmp[9:]
            for s in samples:
                het[s] = {'Hom': 0, 'Het': 0}
        else:
            tmp = line.strip().split()
            genotypes = tmp[9:]
            nsites += 1
            for sam, geno in zip(samples, genotypes):
                calls = geno.split(':')[0]
                #   If the call is missing, then skip it
                if '.' in calls:
                    continue
                elif calls == '0/0' or calls == '1/1':
                    het[sam]['Hom'] += 1
                elif calls == '0/1':
                    het[sam]['Het'] += 1
                else:
                    continue

#   Calculate the heterozygosity
for sample in het.keys():
    ho = float(het[sample]['Het'])/nsites
    print sample, ho
