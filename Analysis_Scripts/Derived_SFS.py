#!/usr/bin/env python

#   Builds the derived SFS, given a VCF of SNP calls and the ancestral state
#   table produced by Get_Ancestral_State.py

import sys

#   the files we operate on
anc = sys.argv[1]
snps = sys.argv[2]

#   a few boolean values that shape how we build the sfs
ignore_lowcoverage = True
ignore_lowqual = False
ignore_hets = True

#   a function to filter our SNPs based on the above flags
def apply_filters(flag, depth, qual, het):
    if depth and flag == 'low_depth':
        return False
    elif qual and flag == 'low_conf':
        return False
    elif het and flag == 'trans_specific':
        return False
    else:
        return True


#   parse the ancestral table into a dictionary
ancestral_states = {}
with open(anc, 'r') as f:
    for line in f:
        tmp = line.strip().split('\t')
        if tmp[0] not in ancestral_states:
            #   check the conditions that we set above
            if apply_filters(tmp[2], ignore_lowcoverage, ignore_lowqual, ignore_hets):
                ancestral_states[tmp[0]] = tmp[1]
            else:
                continue

#   start stepping through the vcf and count up the derived alleles
with open(snps, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            tmp = line.strip().split('\t')
            #   build the position lookup like we have for the ancestral table
            pos = tmp[0] + ':' + tmp[1]
            ref = tmp[3]
            alt = tmp[4]
            #   if the position of the current SNP doesn't have ancestral 
            #   state information, then we skip it
            if pos not in ancestral_states:
                continue
            else:
                #   we get which allele is derived
                #   we use 0 and 1, since that's how they are coded in the VCF
                if ancestral_states[pos] == ref:
                    ancestral = '0'
                    derived = '1'
                elif ancestral_states[pos] == alt:
                    ancestral = '1'
                    derived = '0'
                else:
                    #   if vulgare segregates for a different base, then we 
                    #   have an infinite sites violation, and we skip it
                    continue
                #   Get the genotypes 
                #   the sample information starts on the 10th item
                genotypes = tmp[9:]
                #   we will split the genotype into a list of calls
                calls = []
                for g in genotypes:
                    info = g.split(':')
                    #   we skip over missing genotypes
                    if '.' in info[0]:
                        continue
                    else:
                        #   We split on /, since that's what spearates the
                        #   alleles in the calls
                        calls = calls + info[0].split('/')
                #   Then count up the derived alleles
                derived_count = calls.count(derived)
                #   and print it out
                print pos + '\t' + str(derived_count)
