#!/usr/bin/env python

#   A script to count the number of per-sample deleterious mutations

#   This builds the list of deleterious mutaions
deleterious = []
with open('Deleterious_SIFT.txt', 'r') as f:
    for line in f:
        s = line.strip().split('_')[1]
        deleterious.append(s)

#   This is for barley
samples = ['Barke', 'Bowman', 'Harrington', 'Haruna_Nijo', 'Igri', 'Kindred', 'Morex', 'Steptoe', 'OUH602']
#   This is for soy
#samples = ['Archer', 'Minsoy', 'Noir1', 'Wm82_SGC', 'Williams', 'M92-220', 'IA3023']

for s in samples:
    #   This is the output filename
    filename = s + '_Deleterious.txt'
    #   The BED file with the SNPs per-sample
    snpfilename = s + '_SNPs.bed'
    handle = open(filename, 'w')
    with open(snpfilename, 'r') as f:
        for line in f:
            #   The last column is the ID of the SNP
            snpnumber = line.strip().split('\t')[-1]
            #   We do it this way since this is the only way we can be guaranteed
            #   that we get the correct number
            extension = snpnumber.split('_')[-1]
            if extension in deleterious:
                handle.write(snpnumber + '\n')
            else:
                continue
    handle.close()
