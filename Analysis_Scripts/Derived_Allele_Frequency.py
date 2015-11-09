#!/usr/bin/env python
"""Takes a VCF, GFF, and multiple sequence alignment and uses a "crude" method
to infer the ancestral base at SNP positions listed in the VCF. Uses the
majority state state in the alignment to estimate the ancestral state. Due to
the nature of the multiple sequence alignemnts, this script is designed to only
work on coding variants, but it could be modified to infer noncoding variants
as well.

Arguments: VCF, GFF, MSA directory"""

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio import AlignIO
import gff_parse
import sys
import os
from collections import Counter


class SNP(object):
    """A SNP object. This is just a container for SNP data."""
    def __init__(self):
        self.contig = 'NA'
        self.pos = 'NA'
        self.DAF = 'NA'
        self.ref = 'N'
        self.alt = 'N'
        self.anc = 'N'
        self.snpid = 'NA'
        self.derived_samples = []
        self.genotypes = []

    def get_DAF(self):
        """Calculate the derived allele frequency."""
        #   Figure out which state is the derived allele
        if self.anc == self.ref:
            derived = '1'
        elif self.anc == self.alt:
            derived = '0'
        else:
            derived = 'NA'
        if derived == 'NA':
            self.DAF = 'NA'
        else:
            #   Then, how many derived alleles do we observe?
            total = 0
            der = 0
            for i in self.genotypes:
                if derived in i:
                    der += i.count(derived)
                if not '.' in i:
                    total += 2
            self.DAF = float(der)/total

    def get_derived_samples(self, samples):
        """Given a list of sample names, which samples contain a derived allele
        """
        if self.anc == self.ref:
            derived = '1'
        elif self.anc == self.alt:
            derived = '0'
        else:
            derived = 'NA'
        if derived == 'NA':
            self.derived_samples = ['-']
        else:
            for sam, geno in zip(samples, self.genotypes):
                if derived in geno:
                    self.derived_samples.append(sam)

def read_vcf(vcf_file):
    """Read the VCF and save the SNP position, the reference base, and the
    genotypes."""
    snps = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                tmp = line.strip().split('\t')
                samples = tmp[9:]
            else:
                tmp = line.strip().split('\t')
                #   Start a new SNP object
                curr_snp = SNP()
                #   Extract the position, ID, reference base, and alternate
                #   base for the variant
                curr_snp.contig = tmp[0]
                curr_snp.pos = int(tmp[1])
                curr_snp.snpid = tmp[2]
                curr_snp.ref = tmp[3]
                curr_snp.alt = tmp[4]
                sample_info = tmp[9:]
                for s in sample_info:
                    #   The sample information fields have the form
                    #       N/N:...:...
                    #   We are only interested in the genotype calls, whis is
                    #   the first element before the comma. Split again on
                    #   forward slash to get the individual alleles
                    curr_snp.genotypes.append(''.join(s.split(':')[0].split('/')))
                #   Tacke it onto the list
                snps.append(curr_snp)
    #   And then return it with the samples
    return (samples, snps)


def read_gff(gff_file):
    """Reads and parses the GFF."""
    parsed_gff = gff_parse.GFFHandler()
    parsed_gff.gff_parse(gff_file)
    return parsed_gff


def build_cds_sequences(cds_chunks):
    """Takes a list of CDS annotations and sticks them together."""
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
    #   Check the strand, and sort accordingly. For forward strand features,
    #   we sort from low to high. For reverse strand features, we sort from
    #   high to low.
    if strand == 1:
        feature_parts.sort(key=lambda s: s.start)
    elif strand == -1:
        feature_parts.sort(key=lambda s: s.start, reverse=True)
    #   And then we concatenate them all together
    joint_feature = sum(feature_parts)
    #   Return it!
    return joint_feature


def get_snp_txid(snp, gff_data):
    """Uses the position information and the parsed GFF data to extract the
    transcript ID and CDS position of a given SNP."""
    overlapping_cds = gff_data.overlapping_feature(
        chrom=snp.contig,
        pos=snp.pos,
        feat_type='CDS')
    #   Take the first CDS feature that is overlapping
    if overlapping_cds:
        feat = overlapping_cds[0]
        #   What transcript is it in?
        txid = gff_data.get_parents(feat.ID)
        txid = txid[0].ID
        #   Get the full CDS
        full_cds = gff_data.get_siblings(
            feat.ID,
            feat_type='CDS')
        full_cds = build_cds_sequences(full_cds)
        #   Get the coordinate of the SNP in the CDS
        #   It should be inside the CDS, so we will enumerate from the CDS
        #   start to the CDS end.
        coding = [p for p in range(full_cds.start, full_cds.end+1) if p in full_cds]
        cds_pos = coding.index(int(snp.pos))
        strand = feat.strand
    else:
        txid = None
        cds_pos = None
        strand = None
    return (txid, cds_pos, strand)


def get_majority(msa, seqid, pos):
    """Gets the majority state at an alignment position. Takes a
    MulitpleSequenceAlignemnt object, a sequence ID to search on, and an
    integer position."""
    #   Find which alignment column we want, but keep in mind that we need to
    #   reference the correct sequence and skip gaps. We start counting from
    #   1 because the 'pos' argument is 1-based.
    for index, rec in enumerate(msa):
        if rec.id == seqid:
            ref_seq = index
            break
    key = 1
    for index, base in enumerate(str(msa[ref_seq].seq)):
        if key == pos:
            correct_col = index
            break
        elif base == '-':
            continue
        else:
            key += 1
    #   Then, slice out the correct column
    column = msa[:, correct_col]
    #   Then, we have to remove the query sequence, since we can't use the
    #   sequence itself to infer its own ancestral state
    col = [c for i, c in enumerate(column) if i != ref_seq]
    #   Turn it into a Counter object, which stores counts as a dictionary
    states = Counter(col)
    #   Remove gaps, since a gap can't be ancestral
    del states['-']
    #   Then, get the most frequent base. The most_frequent() method will
    #   give the most common element. In the case of a tie, it will return an
    #   arbitrary one
    if len(states) < 1:
        return('N', 0)
    majority_state = states.most_common(1)[0][0]
    num_seqs = sum(states.values())
    return (majority_state, num_seqs)


def check_env():
    pass


def main():
    gff = read_gff(sys.argv[2])
    sample_names, SNPs = read_vcf(sys.argv[1])
    aln_dir = sys.argv[3]
    for x in SNPs:
        t_id, cds_nuc_pos, strand = get_snp_txid(x, gff)
        if t_id:
            #   Then, we have to build the filename for the multiple sequence
            #   alignment. Replace dots with underscores.
            aln_name = t_id.replace('.', '_')
            #   Read the alignment as a MultipleSequenceAlignemnt object
            aln = os.path.normpath(os.path.join(aln_dir, aln_name) + '.fasta')
            if not os.path.isfile(aln):
                continue
            aln = AlignIO.read(aln, 'fasta')
            #   And get the majority state
            maj, num_spec = get_majority(aln, aln_name, cds_nuc_pos)
            #   If the CDS is on the reverse strand, we have to complement the
            #   ancestral state
            if strand == '-':
                maj_correct = str(Seq(maj).complement())
            else:
                maj_correct = maj
            #   We set the ancestral frequency to be the majority
            x.anc = maj_correct
            #   Then get the dericed allele frequency
            x.get_DAF()
            #   And then the derived samples
            x.get_derived_samples(sample_names)
        else:
            x.anc = 'N'
            x.DAF = 'NA'
            x.derived_samples = ['-']
        print '\t'.join(
            [
            x.snpid,
            x.anc,
            str(x.DAF),
            ','.join(x.derived_samples)]
            )
    return


main()
