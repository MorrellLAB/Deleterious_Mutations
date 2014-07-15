#!/usr/bin/env python

#   A script to separate the sequences from phytozome clusters into
#   multi-sample FASTA files by cluster ID.

import sys

#   A string to hang on to the cluster IDs
old_cluster_id = ''
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        #   Get the ID of the cluster from the sequence name
        if index == 0:
            #   The cluster ID is the second item in the list
            new_cluster_id = line.split('|')[1]
            handle = open(new_cluster_id + '.fasta', 'a')
            handle.write(line)
            old_cluster_id = new_cluster_id
        elif line.startswith('>'):
            new_cluster_id = line.split('|')[1]
            if new_cluster_id == old_cluster_id:
                handle.write(line)
            else:
                handle.close()
                handle = open(new_cluster_id + '.fasta', 'a')
                handle.write(line)
                old_cluster_id = new_cluster_id
        else:
            handle.write(line)
    handle.close()
