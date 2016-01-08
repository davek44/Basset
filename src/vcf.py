#!/usr/bin/env python
from optparse import OptionParser
import sys

import numpy as np
import pandas as pd
import pysam

from dna_io import dna_one_hot

################################################################################
# vcf.py
#
# Methods and classes to support .vcf SNP analysis.
################################################################################

def cap_allele(allele, cap=5):
    ''' Cap the length of an allele in the figures '''
    if len(allele) > cap:
        allele = allele[:cap] + '*'
    return allele


def snps_seq1(snps, genome_fasta, seq_len):
    ''' Produce an array of one hot coded sequences for a list of SNPs.

    Attrs:
        snps [SNP] : list of SNPs
        genome_fasta (str) : genome FASTA file
        seq_len (int) : sequence length to code

    Return:
        seq_vecs (array) : one hot coded sequences surrounding the SNPs
        seqs [str] : raw sequences
        seq_headers [str] : headers for sequences
    '''
    left_len = seq_len/2 - 1
    right_len = seq_len/2

    # open genome FASTA
    genome = pysam.Fastafile(genome_fasta)

    # initialize one hot coded vector list
    seq_vecs_list = []

    # save sequence strings, too
    seqs = []

    # name sequences
    seq_headers = []

    for snp in snps:
        # specify positions in GFF-style 1-based
        seq_start = snp.pos - left_len
        seq_end = snp.pos + right_len + len(snp.ref_allele) - snp.longest_alt()

        # extract sequence as BED style
        seq = genome.fetch(snp.chrom, seq_start-1, seq_end).upper()

        # verify that ref allele matches ref sequence
        seq_ref = seq[left_len:left_len+len(snp.ref_allele)]
        if seq_ref != snp.ref_allele:
            print >> sys.stderr, 'WARNING: skipping %s because reference allele does not match reference genome: %s vs %s' % (snp.rsid, snp.ref_allele, seq_ref)
            continue

        # one hot code ref allele
        seq_vecs_list.append(dna_one_hot(seq[:seq_len], seq_len))
        seqs.append(seq[:seq_len])

        # name ref allele
        seq_headers.append('%s_%s' % (snp.rsid, cap_allele(snp.ref_allele)))

        for alt_al in snp.alt_alleles:
            # remove ref allele and include alt allele
            seq_alt = seq[:left_len] + alt_al + seq[left_len+len(snp.ref_allele):]

            # one hot code
            seq_vecs_list.append(dna_one_hot(seq_alt, seq_len))
            seqs.append(seq_alt)

            # name
            seq_headers.append('%s_%s' % (snp.rsid, cap_allele(alt_al)))

    # stack
    seq_vecs = np.vstack(seq_vecs_list)

    return seq_vecs, seqs, seq_headers


def vcf_snps(vcf_file, index_snp=False, score=False):
    ''' Load SNPs from a VCF file '''
    vcf_in = open(vcf_file)

    # read through header
    line = vcf_in.readline()
    while line[0] == '#':
        line = vcf_in.readline()

    # read in SNPs
    snps = []
    while line:
        snps.append(SNP(line, index_snp, score))
        line = vcf_in.readline()

    return snps


class SNP:
    ''' SNP

    Represent SNPs read in from a VCF file

    Attributes:
        vcf_line (str)
    '''
    def __init__(self, vcf_line, index_snp=False, score=False):
        a = vcf_line.split()
        if a[0].startswith('chr'):
            self.chrom = a[0]
        else:
            self.chrom = 'chr%s' % a[0]
        self.pos = int(a[1])
        self.rsid = a[2]
        self.ref_allele = a[3]
        self.alt_alleles = a[4].split(',')

        self.index_snp = '.'
        if index_snp:
            self.index_snp = a[5]

        self.score = None
        if score:
            self.score = float(a[6])


    def get_alleles(self):
        ''' Return a list of all alleles '''
        alleles = [self.ref_allele] + self.alt_alleles
        return alleles

    def longest_alt(self):
        ''' Return the longest alt allele. '''
        return max([len(al) for al in self.alt_alleles])
