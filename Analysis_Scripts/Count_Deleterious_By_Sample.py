#!/usr/bin/env python
"""Parses the prediction table and counts how many are deleterious by the three
different prediction methods."""

import sys

min_lrt_seq = 10
#   We predicted on 59,865 codons, so that becomes the basis of our multiple
#   testing correction
lrt_sig = 0.05/59865
#   For soy, we counted 64,087 codons
#lrt_sig = 0.05/64087

#   Start dictionaries to hold SNP IDs of deleterious variants carried by each
#   sample.
sift = {}
pph = {}
lrt = {}
intersect = {}
nonsyn = {}


#   Create a function to assess whether or not the LRT prediction turns out to
#   be deleterious or not.
def lrt_del(eff_line, p=lrt_sig, n=min_lrt_seq):
    """Returns True for deleterious predictions and False otherwise. The
    following criteria will return a deleterious prediction:
        - p-value lower than Bonferroni corrected threshold
        - at least 10 sequences represented
        - Masked constraint < 1
        - Rn == 1 (Ref. allele present only in query species)
        - An == 0 (Alt. allele not present in alignment)
    The last two effectively filter for the query species having a truly
    derived allele."""
    fields = eff_line.strip().split()
    try:
        #   For barley the columns are
        #       constraint: 23
        #       p-value: 24
        #       seqcount: 25
        #       rn: 2nd from last
        #       an: last
        #   For soy:
        #       constraint: 20
        #       p-value: 21
        #       seqcount: 22
        #       rn: 2nd from last
        #       an: last
        #constraint = float(fields[22])
        #pval = float(fields[23])
        #nseq = int(fields[24])
        constraint = float(fields[19])
        pval = float(fields[20])
        nseq = int(fields[21])
        rn = int(fields[-2])
        an = int(fields[-1])
    except ValueError:
        return False
    #   Check all the conditions!
    if (pval < p) and (nseq >= n) and (constraint < 1) and (rn == 1 or an == 0):
        return True
    else:
        return False

handle = open(sys.argv[2], 'w')

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
            #   SIFT prediction is the 18th column for barley
            #      15th for soy
            #sift_pred = tmp[17]
            sift_pred = tmp[14]
            #   PPH is the 19th for barley, 16th for soy
            #pph_pred = tmp[18]
            pph_pred = tmp[15]
            #   Start counting the joint number of nonsynonymous SNPs
            if 'Joint' in nonsyn:
                nonsyn['Joint'].append(snpid)
            else:
                nonsyn['Joint'] = [snpid]
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
                if 'Joint' in sift:
                    sift['Joint'].append(snpid)
                else:
                    sift['Joint'] = [snpid]
                for s in samples:
                    if s == '-' or s == 'NA':
                        continue
                    if s in sift:
                        sift[s].append(snpid)
                    else:
                        sift[s] = [snpid]
            if pph_pred == 'deleterious' or sift_pred == 'NONSENSE':
                if 'Joint' in pph:
                    pph['Joint'].append(snpid)
                else:
                    pph['Joint'] = [snpid]
                for s in samples:
                    if s == '-' or s == 'NA':
                        continue
                    if s in pph:
                        pph[s].append(snpid)
                    else:
                        pph[s] = [snpid]
            if lrt_del(line) or sift_pred == 'NONSENSE':
                if 'Joint' in lrt:
                    lrt['Joint'].append(snpid)
                else:
                    lrt['Joint'] = [snpid]
                for s in samples:
                    if s == '-' or s == 'NA':
                        continue
                    if s in lrt:
                        lrt[s].append(snpid)
                    else:
                        lrt[s] = [snpid]
            #   If it's deleterious by all three methods, then we write it into a separate file
            if (sift_pred == 'DELETERIOUS' and pph_pred == 'deleterious' and lrt_del(line)) or sift_pred == 'NONSENSE':
                if 'Joint' in intersect:
                    intersect['Joint'].append(snpid)
                else:
                    intersect['Joint'] = [snpid]
                for s in samples:
                    if s == '-' or s == 'NA':
                        continue
                    if s in intersect:
                        intersect[s].append(snpid)
                    else:
                        intersect[s] = [snpid]
                handle.write(line)
                handle.flush()

handle.close()

#   Print out a nice table
samplenames = sorted(sift.keys())


print "Sample\tSIFT\tPPH\tLRT\tIntersect"
for s in samplenames:
    sift_snps = set(sift[s])
    pph_snps = set(pph[s])
    lrt_snps = set(lrt[s])
    nonsyn_snps = set(nonsyn[s])
    intersection = set(intersect[s])
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
