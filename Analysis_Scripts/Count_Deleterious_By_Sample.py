#!/usr/bin/env python
"""Parses the prediction table and counts how many are deleterious by the three
different prediction methods."""

import sys

min_lrt_seq = 10
#   We predicted on 59,277 codons, so that becomes the basis of our multiple
#   testing correction
lrt_sig = 0.05/59277

#   Start dictionaries to hold SNP IDs of deleterious variants carried by each
#   sample.
sift = {}
pph = {}
lrt = {}
nonsyn = {}

with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split('\t')
            #   We only want to count nonsynonymous SNPs
            if tmp[7] != 'No':
                continue
            #   First save the SNP ID
            snpid = tmp[0]
            #   SIFT prediction is the 18th column
            sift_pred = tmp[17]
            #   PPH is the 19th
            pph_pred = tmp[18]
            #   LRT is the 24th and 25th
            try:
                lrt_pval = float(tmp[23])
            except ValueError:
                lrt_pval = 1
            try:
                lrt_seqnum = int(tmp[24])
            except ValueError:
                lrt_seqnum = 0
            #   Get the sample names, 14th column
            samples = tmp[13].split(',')
            #   Count up nonsynonymous SNPs by sample. These will be derived
            #   nonsynonymous variants
            for s in samples:
                if s == '-' or s == 'NA':
                    continue
                if s in nonsyn:
                    nonsyn[s].append(snpid)
                else:
                    nonsyn[s] = [snpid]
            #   Start counting up the samples
            if sift_pred == 'DELETERIOUS' or sift_pred == 'NONSENSE':
                for s in samples:
                    if s == '-' or s == 'NA':
                        continue
                    if s in sift:
                        sift[s].append(snpid)
                    else:
                        sift[s] = [snpid]
            if pph_pred == 'deleterious' or sift_pred == 'NONSENSE':
                for s in samples:
                    if s == '-' or s == 'NA':
                        continue
                    if s in pph:
                        pph[s].append(snpid)
                    else:
                        pph[s] = [snpid]
            if (lrt_pval <= lrt_sig and lrt_seqnum >= min_lrt_seq) or sift_pred == 'NONSENSE':
                for s in samples:
                    if s == '-' or s == 'NA':
                        continue
                    if s in lrt:
                        lrt[s].append(snpid)
                    else:
                        lrt[s] = [snpid]

#   Create an entry for joint samples
sift['Joint'] = []
for s in sift.iteritems():
    sift['Joint'] += s[1]
pph['Joint'] = []
for s in pph.iteritems():
    pph['Joint'] += s[1]
lrt['Joint'] = []
for s in lrt.iteritems():
    lrt['Joint'] += s[1]
nonsyn['Joint'] = []
for s in nonsyn.iteritems():
    nonsyn['Joint'] += s[1]

#   Print out a nice table
samplenames = sorted(sift.keys())


print "Sample\tSIFT\tPPH\tLRT\tIntersect"
for s in samplenames:
    sift_snps = set(sift[s])
    pph_snps = set(pph[s])
    lrt_snps = set(lrt[s])
    nonsyn_snps = set(nonsyn[s])
    intersection = sift_snps & pph_snps & lrt_snps
    #   Get the lengths
    lsift = len(sift_snps)
    lpph = len(pph_snps)
    llrt = len(lrt_snps)
    lint = len(intersection)
    lnon = len(nonsyn[s])
    toprint = '\t'.join([
        s,
        str(lsift) + ' (' + str(round(lsift/float(lnon), 3)) + ')',
        str(lpph) + ' (' + str(round(lpph/float(lnon), 3)) + ')',
        str(llrt) + ' (' + str(round(llrt/float(lnon), 3)) + ')',
        str(lint) + ' (' + str(round(lint/float(lnon), 3)) + ')'
        ])
    print toprint
