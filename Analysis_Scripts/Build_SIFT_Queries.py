#!/usr/bin/env python

#   A script to build the protein sequences that SIFT needs to assess
#   whether a SNP is deleterious or not
#   takes four arguments:
#       1) The modified GFF that contains the SNP names as the final column
#       2) The VCF file that describes the SNPs
#       3) The annotation output from mRNA_SNPs.py
#       4) The reference FASTA

import sys
import math
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna


sys.stderr.write('Opening '+sys.argv[4]+' for reading...')
#   Read in the reference
handle = open(sys.argv[4], 'r')
ref_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
handle.close()
sys.stderr.write(' Done!\n')

#   Read in the VCF file
sys.stderr.write('Opening ' + sys.argv[2] + ' for reading...')
bed_snps = {}
with open(sys.argv[2], 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            tmp = line.strip().split('\t')
            #   The third column in a VCF is the SNP ID
            bed_snps[tmp[2]] = (tmp[0], tmp[1])
#   Read in the annotation table
#   We only care about this for the amino acid states
gff_annotations = {}
with open(sys.argv[3], 'r') as f:
    for line in f:
        tmp = line.strip().split('\t')
        gff_annotations[tmp[0]] = (tmp[2], tmp[3])

sys.stderr.write('Opening '+sys.argv[1]+' for reading...')
#   Open the GFF describing just CDS with SNPs inside them
#   Extract the CDS as a feature object, so we can get the sequence later
cds_features = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        tmp = line.strip().split('\t')
        chrom = tmp[0]
        snp_id = tmp[-1]
        start = int(tmp[3]) - 1
        end = int(tmp[4]) + 1
        #   There is a 'phase' field for CDS annotations
        #   This describes how many bases must be subtracted from the start
        #   of the feature to reach the beginning of the next codon
        phase = int(tmp[7])
        direction = tmp[6]
        if direction == '-':
            #   If the direction is -, then we take the phase off the end
            strand = -1
            end = end - phase
        else:
            #   If the direction is +, then we take the phase off the front
            #   But 'subtracting' in this direction is adding
            strand = +1
            start = start + phase
        remainder = (end - start) % 3
        #   Now slice off the ends of the incomplete codons
        if direction == '-':
            #   the end has already been altered, so do the start
            start = start + remainder - 1
            #   The -1 strand complicates things with the phase
            #   we have to take one off the end, to correct for python counting
            end = end - 1
        else:
            end = end - remainder
        cds_features[snp_id] = ((chrom, SeqFeature(FeatureLocation(start, end, strand=strand), type='CDS')))
sys.stderr.write(' Done!\n')

for each in bed_snps:
    #   We open handles to the files so we can write the sequences and the
    #   substitutions
    seq_handle = open(each + '_Seq.fasta', 'w')
    sub_handle = open(each + '_Sub.txt', 'w')
    #   The SNP IDs and the mRNA IDs have been associated by zipping
    #   We can get the index of the SNP, and thus get the ID of the mRNA
    #   The tuple looks like
    cds_annotation = cds_features[each][1]
    cds_reference_positions = range(cds_annotation.location.start, cds_annotation.location.end)
    #   Extract the sequence of the mRNA
    #   Biopython provides an extract() method that gets the sequence
    #   of a feature out of its parent
    #   extract() handles strandedness, and will RC something that is on the complement strand
    #   the second element in the tuple is the chromosome
    cds_sequence = Seq(str(cds_annotation.extract(ref_dict[bed_snps[each][0]]).seq))
    #   Translate the sequence into protein
    #   since extract() pulled the right direction and strand, this should
    #   be good enough to get a translation
    protein_sequence = SeqRecord(cds_sequence.translate(), id=each)
    #   Save the start of the CDS feature
    #   These should be in order
    cds_start = cds_annotation.location.start
#        #   Save the start of the mRNA feature
#        mrna_start = mrna_annotation.location.start
    #   Step through the mRNA sequence, saving positions that are in the CDS
#    cds_rel_positions = []
#    for index, base in enumerate(cds_sequence):
#        if index + cds_start in cds_annotation:
#            cds_rel_positions.append(index)
#        #   Get the position of the SNP in the CDS
    snp_cds_position = cds_reference_positions.index(int(bed_snps[each][1]))
    if cds_annotation.strand == 1:
        codon_number = str(int(math.ceil(snp_cds_position / 3)) + 1)
    else:
        codon_number = str(len(protein_sequence) - int(math.ceil(snp_cds_position / 3)))
    states = gff_annotations[each]
    #   Write it out
    SeqIO.write(protein_sequence, seq_handle, 'fasta')
    sub_handle.write(states[0] + codon_number + states[1])
    seq_handle.close()
    sub_handle.close()
