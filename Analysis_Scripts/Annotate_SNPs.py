#!/usr/bin/env python

#   Takes thesearguments
#       1:  BED describing the SNPS
#           Has the form
#           CHR [tab] START [tab] END [tab] ALT_STATE [tab] SNP_ID
#       2:  Reference FASTA
#       3:  CDS with SNPs in them

import sys
#   Load up all our required biopython modules
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna

sys.stderr.write('Opening '+sys.argv[2]+' for reading...')
#   Read in the reference
handle = open(sys.argv[2], 'r')
ref_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
handle.close()
sys.stderr.write(' Done!\n')

sys.stderr.write('Opening '+sys.argv[1]+' for reading...')
#   Read the BED file with the SNPs, save the chromosome, the start, the state
#   and the SNP ID
bed_snps =[]
with open(sys.argv[1], 'r') as f:
    for line in f:
        tmp = line.strip().split('\t')
        snp_id = tmp[4]
        #   Have to cast to integer to do math 
        #   (and checking for set membership in other sets of integers)
        start = int(tmp[1])
        chrom = tmp[0]
        #   Can be a comma-delimited list
        alt = tmp[3].split(',')
        bed_snps.append((snp_id, chrom, start, alt[0]))

#   Open the GFF describing just CDS with SNPs inside them
#   Extract the CDS as a feature object, so we can get the sequence later
cds_snp_ids = []
cds_features = []
with open(sys.argv[3], 'r') as f:
    for line in f:
        tmp = line.strip().split('\t')
        chrom = tmp[0]
        snp_id = tmp[3]
        start = int(tmp[1])
        end = int(tmp[2]) + 1
        #   There is a 'phase' field for CDS annotations
        #   This describes how many bases must be subtracted from the start
        #   of the feature to reach the beginning of the next codon
        phase = int(tmp[5])
        direction = tmp[4]
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
            start = start + remainder
        else:
            end = end - remainder
        cds_features.append((chrom, SeqFeature(FeatureLocation(start, end, strand=strand), type='CDS')))
        cds_snp_ids.append(snp_id)
sys.stderr.write(' Done!\n')
combined = zip(cds_snp_ids, cds_features)

#   Iterate through the list of SNPs from the BED
#   Check if the SNP id is in the list of snp_ids from the mRNA GFF
#   if it is:
#       Extract the position of the SNP from the tuple
#       Extract all CDS annotations under that mRNA
#       Check if the position of the SNP lands in a CDS
#       If it does:
#           Extract the mRNA sequence from the reference
#           enumerate() through the mRNA sequence, saving
#               the index of the bases that land in the CDS (have to check index + mRNA.location.start)
#           Get the index of the SNP position in the CDS annotations
#           Get the position modulus 3
#           Calculate the codon positions
#           Extract the codon from the CDS sequence
#           Copy the codon into a new variable
#           Translate the codons
#           Check if the amino acid states are the same
#           If they are:
#               Print the SNP ID, the position (1, 2, 3), amino acids and 'yes'
#           Else:
#               Print the SNP ID, the position (1, 2, 3), amino acids, and 'no'
#       Else:
#           Print the SNP ID, 'noncoding', two dashes, and 'yes'
#   Else:
#       Print the SNP ID, 'noncoding', two dashes, and 'yes'
for each in bed_snps:
    #   The first element is the SNP ID
    #   and this list is taken from the GFF with mRNAs that contain SNPs
    if each[0] in cds_snp_ids:
        #   If it is in the CDS, then it is coding
        coding = True
        #   The third element is the position
        snp_ref_position = each[2]
    else:
        coding = False
    if coding:
        #   The SNP IDs and the mRNA IDs have been associated by zipping
        #   We can get the index of the SNP, and thus get the ID of the mRNA
        #   The tuple looks like
        #   (SNP ID, (chromosome, feature))
        #   It is ugly, but it maintains all the information
        snpid_index = cds_snp_ids.index(each[0])
        cds_annotation = combined[snpid_index][1][1]
        cds_reference_positions = range(cds_annotation.location.start, cds_annotation.location.end)
        #   Extract the sequence of the mRNA
        #   Biopython provides an extract() method that gets the sequence
        #   of a feature out of its parent
        #   extract() handles strandedness, and will RC something that is on the complement strand
        #   the second element in the tuple is the chromosome
        cds_sequence = cds_annotation.extract(ref_dict[each[1]])
        #   Save the start of the CDS feature
        #   These should be in order
        cds_start = cds_annotation.location.start
        #   Step through the mRNA sequence, saving positions that are in the CDS
        cds_rel_positions = []
        for index, base in enumerate(cds_sequence):
            if index + cds_start in cds_annotation:
                cds_rel_positions.append(index)
        #   Sometimes, the SNPs occur in a junction between CDS
        #   if a codon is split between segments of CDS
        #   We hope this doesn't happen too often, so we'll have to check these
        #   manually. We print a message
        if snp_ref_position not in cds_reference_positions:
            silent = '-'
            translations = ['-', '-']
            codons = ['-', '-']
            position = 'CHECK_ME'
        else:
            #   Get the position of the SNP in the CDS
            snp_cds_position = cds_reference_positions.index(snp_ref_position)
            #   The relative position in the codon
            snp_coding_rel_position = snp_cds_position % 3
            #   The codon positions
            codon_lower_bound = snp_cds_position - snp_coding_rel_position
            codon_upper_bound = snp_cds_position + (2 - snp_coding_rel_position) + 1
            codon_positions = cds_rel_positions[codon_lower_bound:codon_upper_bound]
            #   If the mRNA is on the reverse strand, we should reverse the calls
            #   This is only a problem for 3rd and 1st position SNPs
            #   There is probably an elegant way to do this, but I will use if
            #   since it is easier
            #   It pulls the proper codon, but it doesn't report the correct
            #   position
            if cds_annotation.strand == -1:
                #   1 should be 3
                #   But we are still in python counting
                if snp_coding_rel_position == 0:
                    snp_coding_rel_position = 2
                #   and 3 should be 1
                elif snp_coding_rel_position == 2:
                    snp_coding_rel_position = 0
            #   Set the 'position' variable for output
            position = str(snp_coding_rel_position + 1)
            #   The bases that cover the codon
            ref_codon = [ref_dict[each[1]][pos + cds_start - 1] for pos in codon_positions]
            #   If the strand is reverse, then we should RC the codon and the alt state
            #   The alternate state
            alt_state = each[3]
            if cds_annotation.strand == -1:
                ref_codon = list(Seq(''.join(ref_codon)).reverse_complement())
                alt_state = str(Seq(alt_state).complement())
            #   The alternate codon
            alt_codon = ref_codon[:]
            #   Drop the alternate state into it
            alt_codon[snp_coding_rel_position] = alt_state
            #   Put them in a list, and put them into strings
            codons = [''.join(ref_codon), ''.join(alt_codon)]
            #   Translate them with a list comprehension
            translations = [str(Seq(c, unambiguous_dna).translate(table=1)) for c in codons]
            #   Are the translations the same?
            if translations[0] == translations[1]:
                silent = 'yes'
            else:
                silent = 'no'
    #   If we are not coding
    else:
        silent = 'yes'
        translations = ['-', '-']
        codons = ['-', '-']
        position = 'Not_in_CDS'
    #   Print it all out
    information = [each[0], position, translations[0], translations[1], codons[0], codons[1], silent]
    print '\t'.join(information)
#   The last two else: blocks in the above code block are redundant, but I 
#   left them in for readability. They shoulnd't have much of a performance
#   hit at all.
