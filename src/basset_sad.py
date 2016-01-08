#!/usr/bin/env python
from optparse import OptionParser
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import subprocess

from dna_io import dna_one_hot
import vcf
from util.utils import message

################################################################################
# basset_sad.py
#
# Compute SNP Accessibility Difference scores for SNPs in a VCF file.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <model_th> <vcf_file>'
    parser = OptionParser(usage)
    parser.add_option('--cuda', dest='cuda', default=False, action='store_true', help='Predict on the GPU [Default: %default]')
    parser.add_option('-d', dest='model_hdf5_file', default=None, help='Pre-computed model output as HDF5 [Default: %default]')
    parser.add_option('-f', dest='genome_fasta', default='%s/data/genomes/hg19.fa'%os.environ['BASSETDIR'], help='Genome FASTA from which sequences will be drawn [Default: %default]')
    parser.add_option('-i', dest='index_snp', default=False, action='store_true', help='SNPs are labeled with their index SNP as column 6 [Default: %default]')
    parser.add_option('-l', dest='seq_len', type='int', default=600, help='Sequence length provided to the model [Default: %default]')
    parser.add_option('-m', dest='min_limit', default=0.1, type='float', help='Minimum heatmap limit [Default: %default]')
    parser.add_option('-o', dest='out_dir', default='sad', help='Output directory for tables and plots [Default: %default]')
    parser.add_option('-s', dest='score', default=False, action='store_true', help='SNPs are labeled with scores as column 7 [Default: %default]')
    parser.add_option('-t', dest='targets_file', default=None, help='File specifying target indexes and labels in table format')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide Torch model and VCF file')
    else:
        model_th = args[0]
        vcf_file = args[1]

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    #################################################################
    # prep SNP sequences
    #################################################################
    # load SNPs
    snps = vcf.vcf_snps(vcf_file, options.index_snp, options.score)

    if options.model_hdf5_file is None:
        # get one hot coded input sequences
        seq_vecs, seqs, seq_headers = vcf.snps_seq1(snps, options.genome_fasta, options.seq_len)

        # reshape sequences for torch
        seq_vecs = seq_vecs.reshape((seq_vecs.shape[0],4,1,seq_vecs.shape[1]/4))

        # write to HDF5
        h5f = h5py.File('%s/model_in.h5'%options.out_dir, 'w')
        h5f.create_dataset('test_in', data=seq_vecs)
        h5f.close()


    #################################################################
    # predict in Torch
    #################################################################
    if options.model_hdf5_file is None:
        if options.cuda:
            cuda_str = '-cuda'
        else:
            cuda_str = ''

        options.model_hdf5_file = '%s/model_out.txt' % options.out_dir
        cmd = 'basset_predict.lua -norm %s %s %s/model_in.h5 %s' % (cuda_str, model_th, options.out_dir, options.model_hdf5_file)
        if subprocess.call(cmd, shell=True):
            message('Error running basset_predict.lua', 'error')

    # read in predictions
    seq_preds = []
    for line in open(options.model_hdf5_file):
        seq_preds.append(np.array([float(p) for p in line.split()]))
    seq_preds = np.array(seq_preds)


    #################################################################
    # collect and print SADs
    #################################################################
    if options.targets_file is not None:
        target_labels = [line.split()[0] for line in open(options.targets_file)]

    sad_out = open('%s/sad_table.txt' % options.out_dir, 'w')

    header_cols = ('rsid', 'index', 'score', 'ref', 'alt', 'target', 'ref_pred', 'alt_pred', 'sad')
    print >> sad_out, ' '.join(header_cols)

    # hash by index snp
    sad_matrices = {}
    sad_labels = {}
    sad_scores = {}

    pi = 0
    for snp in snps:
        # get reference prediction
        ref_preds = seq_preds[pi,:]
        pi += 1

        for alt_al in snp.alt_alleles:
            # get alternate prediction
            alt_preds = seq_preds[pi,:]
            pi += 1

            # normalize by reference
            alt_sad = alt_preds - ref_preds
            sad_matrices.setdefault(snp.index_snp,[]).append(alt_sad)

            # label as mutation from reference
            alt_label = '%s_%s>%s' % (snp.rsid, vcf.cap_allele(snp.ref_allele), vcf.cap_allele(alt_al))
            sad_labels.setdefault(snp.index_snp,[]).append(alt_label)

            # save scores
            sad_scores.setdefault(snp.index_snp,[]).append(snp.score)

            # print table lines
            for ti in range(len(alt_sad)):
                if options.index_snp and options.score:
                    cols = (snp.rsid, snp.index_snp, snp.score, vcf.cap_allele(snp.ref_allele), vcf.cap_allele(alt_al), target_labels[ti], ref_preds[ti], alt_preds[ti], alt_sad[ti])
                    print >> sad_out, '%-13s %-13s %5.3f %6s %6s %12s %6.4f %6.4f %7.4f' % cols
                elif options.index_snp:
                    cols = (snp.rsid, snp.index_snp, vcf.cap_allele(snp.ref_allele), vcf.cap_allele(alt_al), target_labels[ti], ref_preds[ti], alt_preds[ti], alt_sad[ti])
                    print >> sad_out, '%-13s %-13s %6s %6s %12s %6.4f %6.4f %7.4f' % cols
                elif options.score:
                    cols = (snp.rsid, snp.score, vcf.cap_allele(snp.ref_allele), vcf.cap_allele(alt_al), target_labels[ti], ref_preds[ti], alt_preds[ti], alt_sad[ti])
                    print >> sad_out, '%-13s %5.3f %6s %6s %12s %6.4f %6.4f %7.4f' % cols
                else:
                    cols = (snp.rsid, vcf.cap_allele(snp.ref_allele), vcf.cap_allele(alt_al), target_labels[ti], ref_preds[ti], alt_preds[ti], alt_sad[ti])
                    print >> sad_out, '%-13s %6s %6s %12s %6.4f %6.4f %7.4f' % cols

    sad_out.close()


    #################################################################
    # plot SAD heatmaps
    #################################################################
    for ii in sad_matrices:
        # convert fully to numpy arrays
        sad_matrix = abs(np.array(sad_matrices[ii]))
        print ii, sad_matrix.shape

        if sad_matrix.shape[0] > 1:
            vlim = max(options.min_limit, sad_matrix.max())
            score_mat = np.reshape(np.array(sad_scores[ii]), (-1, 1))

            if options.targets_file is None:
                # plot heatmap
                plt.figure(figsize=(20, 0.5*sad_matrix.shape[0]))

                # lay out scores
                cols = 12
                ax_score = plt.subplot2grid((1,cols), (0,0))
                ax_sad = plt.subplot2grid((1,cols), (0,1), colspan=(cols-1))

                sns.heatmap(score_mat, xticklabels=False, yticklabels=False, vmin=0, vmax=1, cmap='Reds', cbar=False, ax=ax_score)
                sns.heatmap(sad_matrix, xticklabels=False, yticklabels=sad_labels[ii], vmin=0, vmax=vlim, ax=ax_sad)

            else:
                # plot heatmap
                plt.figure(figsize=(20, 0.5 + 0.5*sad_matrix.shape[0]))

                # lay out scores
                cols = 12
                ax_score = plt.subplot2grid((1,cols), (0,0))
                ax_sad = plt.subplot2grid((1,cols), (0,1), colspan=(cols-1))

                sns.heatmap(score_mat, xticklabels=False, yticklabels=False, vmin=0, vmax=1, cmap='Reds', cbar=False, ax=ax_score)
                sns.heatmap(sad_matrix, xticklabels=target_labels, yticklabels=sad_labels[ii], vmin=0, vmax=vlim, ax=ax_sad)

                for tick in ax_sad.get_xticklabels():
                    tick.set_rotation(-45)
                    tick.set_horizontalalignment('left')
                    tick.set_fontsize(5)

            plt.tight_layout()
            if ii == '.':
                out_pdf = '%s/sad_heat.pdf' % options.out_dir
            else:
                out_pdf = '%s/sad_%s_heat.pdf' % (options.out_dir, ii)
            plt.savefig(out_pdf)
            plt.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
