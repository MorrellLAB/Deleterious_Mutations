#!/usr/bin/env python
"""Script to create a VCF positions file from a SAM file. This script is
designed to work for only the Barley 9K SNP genotyping SNPs against the Morex
draft assembly. Takes three
arguments:
    1) SAM file with contextual sequence mapping data
    2) Illumina design sequences
    3) Reference assembly
This script requires Biopython."""

import sys
import re
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError as e:
    print 'Please install Biopython to use this script.'
    exit(1)

class SNP(object):
    """A class to hold VCF information about an individual SNP. Stores the
    following information:
        SNP ID
        Contig
        Reference Position
        Reference Base
        Alternate Base
        Context Sequence Length
        Reverse Complement?
        Context. Seq. Mapping Location
        The query sequence itself (from SAM)
    """
    def __init__(self):
        self.snpid = None
        self.contig = None
        self.position = None
        self.states = None
        self.reference = None
        self.alternate = None
        self.clen = None
        self.rc = None
        self.mappingloc = None

    def parse_illumina(self, snp_id, illumina):
        """Get the contextual sequence length from the Illumina design sequence.
        Uses a regular expression to search for the '[' that encloses the query
        polymorphism. Also extracts the query states from the Illumina probe
        design."""
        self.snpid = snp_id
        #   Use re.I so that we do case-insensitive match.
        offset = re.search('^[ATCGMRWSYKVHDBN]+\[', illumina, re.I)
        #   .group() returns the matching sequence
        clen = len(offset.group(0))
        self.context_len = clen
        #   We will split on the forward slash (/) to separate the two SNP
        #   states. Some are lowercase for some reason, so we have to force
        #   uppercase
        halves = illumina.split('/')
        self.states = [halves[0][-1].upper(), halves[1][0].upper()]

    def get_mapping_data(self, sam_flag, sam_chr, sam_pos, sam_seq):
        """Set the appropriate mapping data flags, given the information in the
        SAM file. We are interested in the bitwise flag, the contig that the
        context sequence maps to, and the mapping position."""
        #   Check if the 5th (16) bit is set. If so, we need to RC
        if int(sam_flag) & 16:
            self.rc = True
        else:
            self.rc = False
        #   It is easy to get the contig and mapping position
        self.contig = sam_chr
        self.mappingloc = int(sam_pos)
        self.qseq = sam_seq

    def calculate_ref_pos(self):
        """Calculate the reference position on the Morex assembly given the
        contextual sequence length, whether or not the sequence was RCed in
        mapping, and the mapping position of the contextual sequence."""
        if self.rc:
            #   If we reverse-complement, then we have to count from the end
            #   and add to our mapping position accordingly
            self.position = self.mappingloc + \
                len(self.qseq) - (self.context_len)
        else:
            self.position = self.mappingloc + self.context_len - 1

    def get_bases(self, refseq):
        """Check the base in the reference sequence to try to identify the
        correct reference and alternate state."""
        #   Get the actual base out of the reference assembly. Remember to
        #   subtract 1 for 0-based indexing
        morex_base = refseq[self.contig][self.position-1]
        self.reference = morex_base
        #   Then check if the SNP state was RCed when it was mapped. If so,
        #   then we need to complement the allele states before we check for
        #   ref or not.
        if self.rc:
            alleles = [str(Seq(a).complement()) for a in self.states]
        else:
            alleles = self.states
        #   Cast the states to a set so we can easily discard one of them
        alleles = set(self.states)
        alleles.discard(morex_base)
        self.alternate = list(alleles)[0]

    def print_vcf(self):
        """Just print the data in the order that the VCF wants."""
        toprint = [
            self.contig,
            str(self.position),
            self.snpid,
            self.reference,
            self.alternate,
            '.',
            '.',
            's'
            ]
        print '\t'.join(toprint)


def read_reference(refseq):
    """Read the reference sequence, and parse it into a dictionary that we can
    pull out individual contigs."""
    rseq = SeqIO.to_dict(refseq, 'fasta')
    return rseq


def main(sam, illumina, refseq):
    """Main function."""
    #   Read the reference into a dictionary
    reference = SeqIO.to_dict(SeqIO.parse(refseq, 'fasta'))
    #   Define the VCF header
    header = '''##fileformat=VCFv4.2
##INFO<ID=s,Number=1,Type=Flag,Description="Variant is calculated from SAM">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'''
    print header
    #   Start an empty dictionary that will carry our SNP data throughout the
    #   SAM and design sequence parsing
    snp_data = {}
    with open(illumina, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            #   Start a new SNP object
            s = SNP()
            s.parse_illumina(tmp[0], tmp[2])
            snp_data[tmp[0]] = s

    #   Then go through the SAM file
    with open(sam, 'r') as f:
        for line in f:
            #   Skip headers
            if line.startswith('@'):
                continue
            else:
                tmp = line.strip().split()
                #   Unpack the important information
                snpname = tmp[0]
                flag = int(tmp[1])
                ctg = tmp[2]
                map_pos = tmp[3]
                cigar = tmp[5]
                seq = tmp[9]
                #   Check if the contig == *
                #   if it is, then the SNP is unmapped.
                if ctg == '*':
                    continue
                #   Then check the CIGAR string. If it is not only match, then
                #   skip it, too. This is an area that can be improved.
                if not re.match('^[1-9]+M$', cigar):
                    continue
                #   For whatever reason, if the SNP in the SAM is not in the
                #   Illumina design sequence dictionary...
                if snpname not in snp_data:
                    continue
                #   Set the SAM mapping information
                snp_data[snpname].get_mapping_data(
                     flag,
                     ctg,
                     map_pos,
                     seq
                     )
                #   Then calculate the position in the reference
                snp_data[snpname].calculate_ref_pos()
                #   And get the correct bases
                snp_data[snpname].get_bases(reference)
                #   And print out the VCF line
                snp_data[snpname].print_vcf()

#   Run the main function
main(sys.argv[1], sys.argv[2], sys.argv[3])
