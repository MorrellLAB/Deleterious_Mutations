#!/usr/bin/env python

#   A script to create a VCF from the SAM file for the barley iSelect 9k
#   SNP genotyping query sequences. This script will only work for this
#   dataset.
#   Takes three arguments:
#       1) The SAM file
#       2) The FASTA containing the query sequences
#       3) The reference assembly
#   Requires Biopython

import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq

#   Define the IUPAC ambiguities
iupac = {
        'M': ['A', 'C'],
        'R': ['A', 'G'],
        'W': ['A', 'T'],
        'S': ['C', 'G'],
        'Y': ['C', 'T'],
        'K': ['G', 'T']
        }

#   Define the VCF header
header = '''##fileformat=VCFv4.1
##INFO=<IS=s,Number=1,Type=Flag,Description="Variant is calculated from SAM">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'''

#   A function to get the reference and alternate states given the SNP data
def get_states(amb, contig, position):
    if position >= len(ref[contig]):
        return ('N', 'N')
    else:
        reference_state = ref[contig][position]
        #   This is ugly, but it works...
        alt_state = [i for i in iupac[amb] if i != reference_state][0]
        return (reference_state, alt_state)

#   A function to get the actual offset of the SNP, rather than rely 
#   on the expected
def get_context(sequence, expected_context):
    if len(sequence) > expected_context:
        if sequence[expected_context] in iupac and expected_context < len(sequence):
            return (sequence[expected_context], expected_context)
        elif sequence[-expected_context] in iupac and expected_context < len(sequence):
            return (sequence[-expected_context], len(sequence) - expected_context)
        else:
            #   A last resort!
            for i, b in enumerate(sequence):
                if b in iupac:
                    return (b, i)
    else:
        for i, b in enumerate(sequence):
            if b in iupac:
                return (b, i)

#   Read in the FASTA file for the SNPs.
snps = SeqIO.to_dict(SeqIO.parse(sys.argv[2], 'fasta'))

#   Read in the reference sequence, and ocnvert it to a dictionary
#   !!! This takes a lot of time and memory !!!
ref = SeqIO.to_dict(SeqIO.parse(sys.argv[3], 'fasta'))

alignment_data = {}
#   Start reading in the file
with open(sys.argv[1], 'r') as f:
    for line in f:
        #   Skip the sequence header lines
        if line.startswith('@'):
            continue
        else:
            tmp = line.strip().split('\t')
            snp_name = tmp[0]
            #   Get the contextual sequence length
            if snp_name.startswith('11_'):
                clen = 120
            elif snp_name.startswith('12_') or snp_name.startswith('3_'):
                clen = 60
            elif snp_name.startswith('SCRI'):
                clen = 60
            flag = int(tmp[1])
            contig = tmp[2]
            start = int(tmp[3])
            cigar = tmp[5]
            fastaseq = snps[snp_name]
            samseq = tmp[9]
            amb, offset = get_context(fastaseq, clen)
            if snp_name in alignment_data:
                snp_name = snp_name + '_2'
            alignment_data[snp_name] = {'flag':flag, 
                                        'contig':contig, 
                                        'start':start, 
                                        'cigar':cigar, 
                                        'fastaseq':fastaseq, 
                                        'samseq':samseq,
                                        'amb':amb, 
                                        'offset':offset, 
                                        'name':snp_name}

print header
#   Start iterating through the alignment data, and getting the reference
#   and alternate states
for snp_name, snp_data in alignment_data.iteritems():
    #   If the contig is *, then the contextual sequence is unmapped
    if snp_data['contig'] == '*':
        continue
    else:
        #   Count up the soft clipped bases. 
        #   We count on both sides, since soft clipping can happen at the front
        #   or the back end of a read
        #   We don't have to worry as much about hard clipping, since the bases
        #   in SEQ will tell us
        frontclip = re.search('^[1-9]+S', snp_data['cigar'])
        endclip = re.search('[1-9]+S$', snp_data['cigar'])
        #   if there is no clipping, then we get a NoneType object
        if frontclip:
            #   .group() returns the matching part of the string
            #   we slice off the last character since that is S, and will not
            #   cast to integer nicely
            frontclip_offset = int(frontclip.group(0)[:-1])
        else:
            #   If it's NoneType, then we set the offset to 0
            frontclip_offset = 0
        #   And the same for clipping on the end of the read
        if endclip:
            endclip_offset = int(endclip.group(0)[:-1])
        else:
            endclip_offset = 0
        #   Check the flag, if it's 16, then we have to flip the positions
        #   And use the clipping on the end of the read instead of the front
        if snp_data['flag'] == 16 or (snp_data['flag'] - 256) == 16:
            pos = snp_data['start'] + (len(snp_data['fastaseq']) - endclip_offset) - snp_data['offset']
            a = str(Seq(snp_data['amb']).complement())
        else:
            pos = snp_data['start'] + snp_data['offset'] - frontclip_offset
            a = snp_data['amb']
        refbase ,altbase = get_states(a, snp_data['contig'], pos)
        #   A quick sanity check
        if refbase in iupac[a] and altbase in iupac[a]:
            print '\t'.join([snp_data['contig'], str(pos+1), snp_data['name'], refbase, altbase, '.', '.', 's'])
