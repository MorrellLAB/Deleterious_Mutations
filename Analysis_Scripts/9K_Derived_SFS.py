#!/usr/bin/env python

#   A script to calculate the derived frequency spectrum for the 9K SNPs
#   in the NSGC genotyping matrix from Ana
#   Takes two arguments
#       1) The inferred ancestral states for the 9K SNPs
#       2) The genoyping matrix, which is exported as csv
#           Keep in mind that this matrix has been converted to the BOPAC
#           naming convention; it's not the raw matrix

import sys

#   Open the file containing the inferred ancestral state and read it into
#   a dictionary
ancestral = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        tmp = line.strip().split('\t')
        ancestral[tmp[0]] = tmp[1]

#   Open the genotype matrix and start calculating the derived allele frequency
with open(sys.argv[2], 'r') as f:
    for index, line in enumerate(f):
        #   Skip the first line, which is a header
        if index == 0:
            continue
        else:
            tmp = line.strip().split(',')
            #   The SNP name is the first column in the table
            name = tmp[0]
            #   Just a sanity check - if the SNP name doesn't have an ancestral
            #   state, we skip it
            if name not in ancestral:
                continue
            else:
                #   Put all the calls together into a string; we can call the
                #   count() method on it later
                calls = ''.join(tmp[1:])
                #   If, for some reason the ancestral base isn't segregating in the NSGC, we skip over it
                if ancestral[name] not in calls and len(set(calls)) > 0:
                    continue
                else:
                    #   Get the state that isn't the ancestral state
                    states = set(calls)
                    #   but remove Ns
                    if 'N' in states:
                        states.remove('N')
                    derived_allele = [s for s in list(states) if s != ancestral[name]]
                    if len(derived_allele) == 0:
                        continue
                    else:
                        derived_allele = derived_allele[0]
                    #   Count up the number of derived alleles
                    #   And get the frequency
                    derived_freq = float(calls.count(derived_allele))/len(calls)
                    #   Print it all out!
                    print name + '\t' + str(derived_freq)
