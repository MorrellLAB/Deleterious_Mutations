#!/usr/bin/env python
#   A script to take the putative ancestral states from a set of species and 
#   build a bunch of multiple-sequence alignment files out of them for
#   inference of maximum-likelihood ancestral state with FastML

import sys
from Bio import SeqIO

#   How many contigs are present?
#n_contigs = 2670738
n_contigs=7

#   Get a list of all the arguments passed
arguments = sys.argv[1:]
#   This is a very specific case, but the file names are the species names
species = [s.replace('.fasta', '') for s in arguments]
#   Next, we open generators to each of the fasta files
#       We just trust that they are actually fasta.
handles = [SeqIO.parse(open(f, 'r'), 'fasta') for f in arguments]
#   We have a list for data we have in "waiting"
#   This is for when the files get out of phase and we have to pull sequences
#   from some files but not others
waiting = [None] * len(handles)

#   What we will do is iterate over every contig, and check if the current one
#   is the correct contig. If it is, we use it and pull the next. If not, then
#   we wait
for contig in xrange(1, n_contigs+1):
    expected_contig = 'morex_contig_' + str(contig)
    sequences = []
    for i, h in enumerate(handles):
        #   First, check if we have any data 'in waiting'
        if waiting[i]:
            #   If there is data, then we check if it is the right data
            if waiting[i].name == expected_contig:
                #   If it's the right data, we pull it and remove the data 
                #   from waiting
                seqdata = str(waiting[i].seq)
                waiting[i] = None
            else:
                seqdata = 'N'
        else:
            #   First, pull the next contig from the file
            try:
                s = h.next()
            except StopIteration:
                print 'Done iterating through ' + species[i]
            #   Now that we have pulled some sequence, we have to check which
            #   contig it came from
            if s.name == expected_contig:
                #   If the contig matches the one we expected, then tack it
                #   onto the list
                seqdata = str(s.seq)
            else:
                #   if it doesn't match, then we stick the data in waiting
                waiting[i] = s
                seqdata = 'N'
        sequences.append(seqdata)
    #   If we have any missing, then we want to replace it with Ns of the
    #   correct length
    aligned_len = max([len(x) for x in sequences])
    #   Check for things that are Ns
    sequences = ['N'*aligned_len if x == 'N' else x for x in sequences]
    #   Next we print out the alignments
    handle = open(expected_contig + '.fasta', 'w')
    for sp, aln in zip(species, sequences):
        handle.write('>' + sp + '\n' + aln + '\n')
    handle.close()
