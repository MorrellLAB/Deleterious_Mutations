#!/usr/bin/env python

#   To take arguments on the command line
import sys
#   To operate on lists of numbers, to return intervals
import itertools
import operator
#   To read sequences
from Bio import SeqIO

#   Get a list of the sequences
full_seqs = list(SeqIO.parse(sys.argv[1], 'fasta'))

#   Strip out just the sequence data
#   and convert to a list
raw_seqs = [list(str(s.seq)) for s in full_seqs]
#   Transpose them, to iterate over alignment columns
columns = zip(*raw_seqs)

indel_positions = []
#   For each column
for index, base in enumerate(columns):
    #   Convert to a set, which gives us how many states in each column
    tmp = set(base)
    #   Then we throw out Ns
    tmp.discard('N')
    #   this is a little ugly, but it works
    #   start a new list for the indel positions
    if '-' in tmp:
        indel_positions.append(index+1)

#   This identifies runs of consecutive integers, and puts them into separate lists
#   for group_key, group in (0, pos1) (1, pos2) (2, pos3) ...
#   group by the value of (pos - index) 
#   if they are the same, then the numbers are consecutive
indels = []
for key, group in itertools.groupby(enumerate(indel_positions), lambda (i,x):i-x):
    #   save the second value, which is the actual group
    indels.append(map(operator.itemgetter(1), group))

#   iterate through the list of lists and print out the ranges
for i in indels:
    if len(i) == 1:
        print str(i[0]) + '\t' + str(i[0])
    else:
        print str(i[0]) + '\t' + str(i[-1])
