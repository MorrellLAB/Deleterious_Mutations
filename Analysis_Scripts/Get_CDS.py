#!/usr/bin/env python
"""Extract the CDS sequence of a transcript give its
    - Stable transcript ID
    - GFF (v3) describing its features on the genome
    - Reference FASTA file
Requires Biopython.
"""

#   Load up the modules
import sys
from Bio import SeqIO
#   Required for handling sequence features and performing extractions
from Bio.SeqFeature import SeqFeature, FeatureLocation
#   My class to parse GFF data
import gff_parse


#   Define some functions to handle our data
def read_gff(gff_file):
    """Reads the GFF file and parses it into an easy-to-access data structure.
    Depends on the gff_parse class that I wrote."""
    parsed_gff = gff_parse.GFFHandler()
    parsed_gff.gff_parse(gff_file)
    return parsed_gff


def read_reference(refseq):
    """Reads the reference sequence into a dictionary and returns it. For
    the barley genome, this will be slow and memory-intensive, but should work
    for more nicely-formatted genomes."""
    handle = open(refseq, 'r')
    #   Start a parser that will iterate through the sequence file
    ref_parser = SeqIO.parse(handle, 'fasta')
    #   And then turn the parser into a sequence dictionary
    ref_dict = SeqIO.to_dict(ref_parser)
    return ref_dict


def build_cds_sequences(cds_chunks):
    """Takes a list of CDS annotations and sticks them together, respecting
    their phase fields."""
    #   Next, we build the compound feature that is the CDS
    feature_parts = []
    for c in cds_chunks:
        #   We convert all the coordinates to 0-based, since they came out of
        #   the GFF as 1-based. SeqFeature.extract() will give an interval
        #   corresponding to [start, end), so we only have to subtract from
        #   the start.
        start = int(c.start) - 1
        end = int(c.end)
        if c.strand == '+':
            strand = 1
        elif c.strand == '-':
            strand = -1
        #   Then, we build a FeatureLocation object out of it, and give it the
        #   proper strandedness
        feature_sub = FeatureLocation(start, end)
        feature_parts.append(feature_sub)
        #   And lastly, we save the contig it's on
        ctg = c.seqid
    #   Check the strand, and sort accordingly. For forward strand features,
    #   we sort from low to high. For reverse strand features, we sort from
    #   high to low.
    if strand == 1:
        feature_parts.sort(key=lambda s: s.start)
    elif strand == -1:
        feature_parts.sort(key=lambda s: s.start, reverse=True)
    #   And then we concatenate them all together
    joint_feature = sum(feature_parts)
    joint_feature = SeqFeature(joint_feature, type='CDS', strand=strand)
    #   Return it!
    return(ctg, joint_feature)


def extract_seq(feature, sequence):
    """Extracts a CDS sequence from the reference assembly, given a feature
    for the CDS and a Seq object for the reference."""
    #   Easy-peasy, just call extract()
    s = feature.extract(sequence)
    return s


def usage():
    """Throw a usage message. This will be printed if the script is called
    without any arguments."""
    message = """
Usage: Get_Barley_CDS_Seqs.py [TranscriptIDs] [GFF] [Reference]

Will create a series of single-sequence FASTA files with the CDS sequences for
each of the transcripts listed in [TranscriptIDs]. This script expects the IDs
of the transcripts to be in a plain text file, one per line. The GFF should be
a GFF v3 file, and the reference should be FASTA.

This script requires Biopython."""
    print message
    exit(1)


#   The main function
def main():
    """Main function."""
    #   Check for improper number of arguments. There should be 3
    if len(sys.argv[1:]) != 3:
        usage()
    #   Define our arguments
    transcript_ids = sys.argv[1]
    gff = sys.argv[2]
    reference = sys.argv[3]
    #   Read the transcript IDs into a list
    tx_ids = []
    with open(transcript_ids, 'r') as f:
        for line in f:
            tx_ids.append(line.strip())
    #   Then, we parse the gff and get it into a nice structure
    gff_data = read_gff(gff)
    ref_dict = read_reference(reference)
    #   Next, we get the CDS features for each of our transcripts of interest
    for t in tx_ids:
        cds = gff_data.get_children(t, feat_type='CDS')
        contig, full_cds = build_cds_sequences(cds)
        cds_sequence = extract_seq(full_cds, ref_dict[contig])
        #   We should set the name of the sequence so it will appear in our
        #   FASTA file. But, we want to replace the periods ('.') with
        #   underscores, since they break other analysis tools.
        seqname = t.replace('.', '_')
        cds_sequence.id = seqname
        cds_sequence.description = seqname
        #   Then, we want to create some FASTA files
        handle = open(seqname + '.fasta', 'w')
        SeqIO.write(cds_sequence, handle, 'fasta')
        handle.flush()
        handle.close()

#   Do work
main()
