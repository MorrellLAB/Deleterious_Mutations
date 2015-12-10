#!/usr/bin/env python
"""Calculate per-sample observered heterozygosity as the proportion of all
polymorphic sites that are heterozygous in a VCF. Applies filters on
missingness, read depth, and allelic balance before counting a variant as a
'true' heterozygous call."""

import sys

#   Our thresholds for filtering false positive heterozygous sites
#       Mininum depth for a single sample. 15 reads minimum, 100 reads max.
mindp = 15
maxdp = 100
#       Like "minor allele frequency" - the minimum amount of deviation from
#       50-50 ref/alt reads. +/- 10% deviation from 50-50
mindev = 0.10
#       The maximum amount of missingness to count a site, across samples.
#       Anything with more than 3 samples missing (~25%) will be filtered
maxmis = 3
#   Empty dictionary to store the numbers of homozygous and heterozygous
#   sites
het = {}

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
            #   These are the field we will use to filter our genomic info:
            #       DP (3rd field) for depth
            #       AD (5th field) for allele balace
            genotypes = tmp[9:]
            #   Check the missing threshold
            miss = 0
            for i in genotypes:
                if i.startswith('./.'):
                    miss += 1
            if miss > maxmis:
                continue
            else:
                for sam, geno in zip(samples, genotypes):
                    gt_metadata = geno.split(':')
                    #   Skip missing sites
                    if './.' in gt_metadata or '.' in gt_metadata:
                        continue
                    calls = gt_metadata[0]
                    dp = gt_metadata[2]
                    ad = gt_metadata[4].split(',')
                    #   If the depth is not sufficient, skip it. Also if the
                    #   depth is TOO high, skip it (potential paralogue)
                    if int(dp) < mindp or int(dp) > maxdp:
                        continue
                    elif calls == '0/0' or calls == '1/1':
                        het[sam]['Hom'] += 1
                    elif calls == '0/1':
                        #   If the allelic balance is too skewed, skip it.
                        #   likely sequencing error or paralogous alignment
                        ref = float(ad[0])
                        alt = float(ad[1])
                        balance = ref/(ref+alt)
                        if abs(0.5-balance) > mindev:
                            continue
                        else:
                            het[sam]['Het'] += 1
                    else:
                        continue

#   Calculate the heterozygosity
for sample in het.keys():
    ho = float(het[sample]['Het'])/(het[sample]['Hom'] + het[sample]['Het'])
    print sample + '\t' + str(ho)
