#!/usr/bin/env python
"""Convert from genome-based coordinates to gene-based coordinates for use with
the LRT_Predict pipeline. This has very sepcific file format requirements:

A well-formated GFF file with MLOC parent ID in the final column, and
A BED file describing variants
"""

import sys
import math

gff = sys.argv[1]
bed = sys.argv[2]

#   Read in the GFF data as a dictionary
gff_data = {}
with open(gff, 'r') as f:
    for line in f:
        tmp = line.strip().split()
        #   Get the MLOC parent ID
        #   This is a little gross because the column is formatted like
        #       ID=blahblah;parent=blahblah
        parent = tmp[-1].split(';')[1].split('=')[1]
        #   And the contig
        contig = tmp[0]
        #   Next, the strand
        strand = tmp[6]
        if parent not in gff_data:
            #   If the parent ID is not in the dictionary yet, we start a new
            #   entry for it
            gff_data[parent] = {
                'contig': contig,
                'strand': strand,
                'positions': []
            }
        #   Then, start building the positions covered by the feature
        #   First, the start and stop.
        start = int(tmp[3])
        stop = int(tmp[4])
        #   And the phase
        phase = int(tmp[7])
        #   Then, adjust the start and stop coordinates by the phase as
        #   necessary. According to the GFF3 documentation, if the strand is +
        #   then we subtract the phase from the start, and if the strand is -,
        #   then we subtract the phase from the stop.
        if strand == '-':
            stop = stop - phase
        elif strand == '+':
            start = start - phase
        #   We then generate a list of positions that are covered by the CDS.
        segment = range(start, stop)
        gff_data[parent]['positions'] += segment

#   Next, we go through the BED file and get the gene ID, residue number and
#   SNP ID for each of the SNPs in our testing dataset.
with open(bed, 'r') as f:
    for line in f:
        tmp = line.strip().split()
        snpid = tmp[-1]
        #   We want the 1-based coordinate for the SNP, so we take the end
        pos = int(tmp[2])
        #   This is a little bit messy, but we need to use the contig ID to
        #   match the SNP with the gene name.
        contig = tmp[0]
        #   Iterate through the GFF dictionary and get the gene ID
        for geneid, gffdata in gff_data.iteritems():
            if gffdata['contig'] == contig:
                #   Get the position in the CDS
                cds_pos = gffdata['positions'].index(pos)
                if gffdata['strand'] == '-':
                    cds_pos = len(gffdata['positions']) - cds_pos
                #   Then, get the residue number
                residue = math.floor(cds_pos / 3) + 1
                #   Now, print it all out
                print '\t'.join([geneid, str(int(cds_pos)), snpid])
