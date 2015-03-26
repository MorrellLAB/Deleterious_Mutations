#!/usr/bin/env python

#   A script to fix the ancestral sequence from FastML output
#   For some reason, when all sequences in the alignment have 'N' FastML will
#   output 'G' as the ancestral base. I do not want this. This script will
#   parse the source alignment and the ancestral sequence, replacing all sites
#   that are 'N' in the alignment with 'N' in the ancestral sequence.
#   Arguments are:
#       1) Alignment of ancestral sequences
#       2) ML Ancestral sequence FASTA file


#   Gotta take arguments
import sys

#   To parse alignments
from Bio import AlignIO
#   and sequences
from Bio import SeqIO
#   These let us modify a sequence in-place
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna

#   What is the output filename?
outfile = sys.argv[2].replace('.fasta', '_Fixed.fasta')

#   Read in the alignment
#   These should be small so we can read them all up at once
anc_align = AlignIO.read(sys.argv[1], 'fasta')
raw_ancestral = SeqIO.read(sys.argv[2], 'fasta')

#   We have to create a MutableSeq object, since normal Seq objects are
#   immutable. generic_dna lets us know that we are working with nucleotide
#   sequences
mut_ancestral = MutableSeq(str(raw_ancestral.seq), generic_dna)

#   Next, we iterate through the alignment, column-by-column, and check if
#   all the sequences are 'N.' If they are, then we change the ML ancestral
#   sequence to 'N' as well
#       There is not direct way to do that, so we index by column instead with
#       range(). We just take the first sequence since they should all be the
#       same length in an alignment anyway
for i in range(0, len(anc_align[0])):
    #   This construct extracts column i from the alignment as a string
    column = anc_align[:, i]
    #   Check the number of 'N' in the column
    if column.count('N') == len(column):
        #   Replace the corresponding position in the ancestral sequence
        mut_ancestra[i] = 'N'
    else:
        continue
