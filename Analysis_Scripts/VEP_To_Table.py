#!/usr/bin/env python
"""Script to use a VCF with ID and the output from the Ensembl Variant Effect
Predictor (VEP) to create a table that is of the same format as the barley SNP
effect predictions. This is only for the soybean samples in Bad Mutations 1."""

import sys


def parse_vcf(fname):
    """Read in the VCF and store it as a big dictionary. We are only interested
    in saving the position and ID information for now."""
    vcf_data = {}
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                snpid = tmp[2]
                chrom = tmp[0]
                pos = tmp[1]
                vcf_data[(chrom, pos)] = snpid
    return vcf_data


def parse_vep(vep_line, vcf_data):
    """Read through the VEP table, and associate the position information with
    SNP ID from the VCF."""
    tmp = vep_line.strip().split()
    #   The 'position' has all the info we need.
    pos = tuple(tmp[1].split(':'))
    if pos in vcf_data:
        snpid = vcf_data[pos]
    else:
        return None
    #   Now that we have the SNP ID, we start compiling the other
    #   info we are interested in. The information we will store,
    #   in order is:
    #       1: SNP ID
    #       2: Chromosome
    #       3: Reference Position
    #       4: Transcript ID
    #       5: Reference Base
    #       6: Alternate Base
    #       7: CDS Position
    #       8: Silent (Yes/No)
    #       9: Codon Position
    #       10: AA1 (Reference AA)
    #       11: AA2 (Alternate AA)
    chrom = pos[0]
    position = pos[1]
    txid = tmp[4]
    refbase = tmp[0].split('/')[0][-1]
    altbase = tmp[0].split('/')[1]
    cds_pos = tmp[9].split('/')[0]
    if 'missense' in tmp[6] or 'stop_gained' in tmp[6]:
        silent = 'No'
    else:
        silent = 'Yes'
    codon = tmp[11].split('/')[0]
    codon_pos = '-'
    for index, c in enumerate(codon):
        if c.isupper():
            codon_pos = str(index+1)
    if tmp[10] == '-':
        aa1 = '-'
        aa2 = '-'
    elif len(tmp[10]) == 1:
        aa1 = tmp[10]
        aa2 = tmp[10]
    else:
        aa1, aa2 = tmp[10].split('/')
    return (
        snpid,
        chrom,
        position,
        txid,
        refbase,
        altbase,
        cds_pos,
        silent,
        codon_pos,
        aa1,
        aa2
        )


def main():
    """Main function."""
    vdata = parse_vcf(sys.argv[1])
    #   Print the header"
    print '\t'.join([
        'SNPID',
        'Chromosome',
        'Position',
        'TranscriptID',
        'RefBase',
        'AltBase',
        'CDS_Pos',
        'Silent',
        'CodonPosition',
        'AA1',
        'AA2'
        ])
    with open(sys.argv[2], 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                vepline = parse_vep(line, vdata)
                if vepline:
                    toprint = '\t'.join(list(vepline))
                    print toprint

main()
