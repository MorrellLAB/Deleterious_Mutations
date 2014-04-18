#!/usr/bin/env python

#   A script to extract the H. bulbosum state from a VCF of H. bulbosum aligned
#   against H. vulgare. Assumes a 1-sample VCF with genotype calls
#   ignores indels
#   Note that since we have a sample size of 2 chromosomes, we will not have a 
#   good estimate of trans-species segregation between bulbosum and vulgare

import sys

#   Print the header
print 'Position\tAncestral\tFlag'
#   open and start parsing the VCF
with open(sys.argv[1], 'r') as f:
    for line in f:
        #   Skip lines that start with #
        if line.startswith('#'):
            continue
        else:
            #   Strip the newline and split the line on tabs
            tmp = line.strip().split('\t')
            #   The sample information is the 10th element of the list
            sample = tmp[9]
            #   We also want the contig, position, ref and alt alleles
            contig = tmp[0]
            pos = tmp[1]
            qual = tmp[5]
            ref = tmp[3]
            alt = tmp[4]
            #   We can't call indels very well
            if len(ref) > 1 or len(alt) > 1:
                continue
            #   The genotype is contained in the first part of the sample
            #   information field
            genotype = sample.split(':')[0]
            #   if we have a missing genotype, we just skip it
            if '.' in genotype:
                continue
            #   we also want the depth
            depth = sample.split(':')[2]
            #   get the calls
            calls = genotype.split('/')
            #   If the call is heterozygous, then we make a note
            if len(set(calls)) > 1:
                flag = 'trans-specific'
            elif qual == '.' or float(qual) < 10.0:
                flag = 'low_conf'
            elif int(depth) < 5:
                flag = 'low_depth'
            else:
                flag = 'inferred_ancestral'
            #   then we use the genotype to get the inferred ancestral state
            if len(set(calls)) > 1:
                ancestral = ref+alt
            elif genotype == '0/0':
                ancestral = ref
            elif genotype == '1/1':
                ancestral = alt
            else:
                ancestral = 'N'
            print contig + ':' + pos + '\t' + ancestral + '\t' + flag
