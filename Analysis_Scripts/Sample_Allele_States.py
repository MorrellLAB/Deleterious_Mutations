#!/usr/bin/env python
"""Script to print out how many samples carry the reference and the alternate
alleles at a locus. Will not count or print missing calls. Will also be
updated in the future to report ancestral and derived."""

import sys

print '\t'.join(['SNP_ID', 'Ref', 'Alt', 'Num_Ref', 'Num_Het', 'Num_Alt', 'Ref_Samples', 'Het_Samples', 'Alt_Samples'])

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            samples = line.strip().split()[9:]
        else:
            tmp = line.strip().split()
            snp_id = tmp[2]
            ref = tmp[3]
            alt = tmp[4]
            genotype_fields = tmp[9:]
            ref_samples = []
            het_samples = []
            alt_samples = []
            for g, s in zip(genotype_fields, samples):
                call = g.split(':')[0]
                states = call.split('/')
                if len(set(states)) == 1:
                    if states[0] == '0':
                        ref_samples.append(s)
                    elif states[0] == '1':
                        alt_samples.append(s)
                elif len(set(states)) == 2:
                    het_samples.append(s)
            num_ref = len(ref_samples)
            num_het = len(het_samples)
            num_alt = len(alt_samples)
            if len(ref_samples) == 0:
                ref_samples = ['-']
            if len(het_samples) == 0:
                het_samples = ['-']
            if len(alt_samples) == 0:
                alt_samples = ['-']
            toprint = '\t'.join([snp_id, ref, alt, str(num_ref), str(num_het), str(num_alt), ','.join(ref_samples), ','.join(het_samples), ','.join(alt_samples)])
            print toprint
