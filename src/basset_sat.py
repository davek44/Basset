#!/usr/bin/env python
from optparse import OptionParser
import copy, os, pdb, random, sys

import h5py
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import seaborn as sns

import dna_io
from seq_logo import seq_logo
import subprocess
from util.utils import message

################################################################################
# basset_sat.py
#
# Perform an in silico saturated mutagenesis of the given test sequences using
# the given model.
################################################################################

################################################################################
# main
############################s####################################################
def main():
    usage = 'usage: %prog [options] <model_file> <input_file>'
    parser = OptionParser(usage)
    parser.add_option('-a', dest='input_activity_file', help='Optional activity table corresponding to an input FASTA file')
    parser.add_option('-d', dest='model_hdf5_file', default=None, help='Pre-computed model output as HDF5 [Default: %default]')
    parser.add_option('-g', dest='gain_height', default=False, action='store_true', help='Nucleotide heights determined by the max of loss and gain [Default: %default]')
    parser.add_option('-m', dest='min_limit', default=0.1, type='float', help='Minimum heatmap limit [Default: %default]')
    parser.add_option('-n', dest='center_nt', default=200, type='int', help='Center nt to mutate and plot in the heat map [Default: %default]')
    parser.add_option('-o', dest='out_dir', default='heat', help='Output directory [Default: %default]')
    parser.add_option('-s', dest='sample', default=None, type='int', help='Sample sequences from the test set [Default:%default]')
    parser.add_option('-t', dest='targets', default='0', help='Comma-separated list of target indexes to plot (or -1 for all) [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide Basset model file and input sequences (as a FASTA file or test data in an HDF file')
    else:
        model_file = args[0]
        input_file = args[1]

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    #################################################################
    # parse input file
    #################################################################
    try:
        # input_file is FASTA

        # load sequences and headers
        seqs = []
        seq_headers = []
        for line in open(input_file):
            if line[0] == '>':
                seq_headers.append(line[1:].rstrip())
                seqs.append('')
            else:
                seqs[-1] += line.rstrip()

        model_input_hdf5 = '%s/model_in.h5'%options.out_dir

        if options.input_activity_file:
            # one hot code
            seqs_1hot, targets = dna_io.load_data_1hot(input_file, options.input_activity_file, mean_norm=False, whiten=False, permute=False, sort=False)

            # read in target names
            target_labels = open(options.input_activity_file).readline().strip().split('\t')

        else:
            # load sequences
            seqs_1hot = dna_io.load_sequences(input_file, permute=False)
            targets = None
            target_labels = None

        # sample
        if options.sample:
            sample_i = np.array(random.sample(xrange(seqs_1hot.shape[0]), options.sample))
            seqs_1hot = seqs_1hot[sample_i]
            seq_headers = seq_headers[sample_i]
            if targets is not None:
                targets = targets[sample_i]

        # reshape sequences for torch
        seqs_1hot = seqs_1hot.reshape((seqs_1hot.shape[0],4,1,seqs_1hot.shape[1]/4))

        # write as test data to a HDF5 file
        h5f = h5py.File(model_input_hdf5, 'w')
        h5f.create_dataset('test_in', data=seqs_1hot)
        h5f.close()

    except (IOError, IndexError):
        # input_file is HDF5

        try:
            model_input_hdf5 = input_file

            # load (sampled) test data from HDF5
            hdf5_in = h5py.File(input_file, 'r')
            seqs_1hot = np.array(hdf5_in['test_in'])
            targets = np.array(hdf5_in['test_out'])
            try: # TEMP
                seq_headers = np.array(hdf5_in['test_headers'])
                target_labels = np.array(hdf5_in['target_labels'])
            except:
                seq_headers = None
                target_labels = None
            hdf5_in.close()

            # sample
            if options.sample:
                sample_i = np.array(random.sample(xrange(seqs_1hot.shape[0]), options.sample))
                seqs_1hot = seqs_1hot[sample_i]
                seq_headers = seq_headers[sample_i]
                targets = targets[sample_i]

                # write sampled data to a new HDF5 file
                model_input_hdf5 = '%s/model_in.h5'%options.out_dir
                h5f = h5py.File(model_input_hdf5, 'w')
                h5f.create_dataset('test_in', data=seqs_1hot)
                h5f.close()

            # convert to ACGT sequences
            seqs = dna_io.vecs2dna(seqs_1hot)


        except IOError:
            parser.error('Could not parse input file as FASTA or HDF5.')


    #################################################################
    # Torch predict modifications
    #################################################################
    if options.model_hdf5_file is None:
        options.model_hdf5_file = '%s/model_out.h5' % options.out_dir
        torch_cmd = 'basset_sat_predict.lua -center_nt %d %s %s %s' % (options.center_nt, model_file, model_input_hdf5, options.model_hdf5_file)
        if subprocess.call(torch_cmd, shell=True):
            message('Error running basset_sat_predict.lua', 'error')

    #################################################################
    # load modification predictions
    #################################################################
    hdf5_in = h5py.File(options.model_hdf5_file, 'r')
    seq_mod_preds = np.array(hdf5_in['seq_mod_preds'])
    hdf5_in.close()

    # trim seqs to match seq_mod_preds length
    seq_len = len(seqs[0])
    delta_start = 0
    delta_len = seq_mod_preds.shape[2]
    if delta_len < seq_len:
        delta_start = (seq_len - delta_len)/2
        for i in range(len(seqs)):
            seqs[i] = seqs[i][delta_start:delta_start+delta_len]

    # decide which cells to plot
    if options.targets == '-1':
        plot_targets = xrange(seq_mod_preds.shape[3])
    else:
        plot_targets = [int(ci) for ci in options.targets.split(',')]


    #################################################################
    # plot
    #################################################################
    table_out = open('%s/table.txt' % options.out_dir, 'w')

    rdbu = sns.color_palette("RdBu_r", 10)

    nts = 'ACGT'
    for si in range(seq_mod_preds.shape[0]):
        try:
            header = seq_headers[si]
        except TypeError:
            header = 'seq%d' % si
        seq = seqs[si]

        # plot some descriptive heatmaps for each individual cell type
        for ci in plot_targets:
            seq_mod_preds_cell = seq_mod_preds[si,:,:,ci]
            real_pred_cell = get_real_pred(seq_mod_preds_cell, seq)

            # compute matrices
            norm_matrix = seq_mod_preds_cell - real_pred_cell
            min_scores = seq_mod_preds_cell.min(axis=0)
            max_scores = seq_mod_preds_cell.max(axis=0)
            minmax_matrix = np.vstack([min_scores - real_pred_cell, max_scores - real_pred_cell])

            # prepare figure
            sns.set(style='white', font_scale=0.5)
            sns.axes_style({'axes.linewidth':1})
            heat_cols = 400
            sad_start = 1
            sad_end = 323
            logo_start = 0
            logo_end = 324
            fig = plt.figure(figsize=(20,3))
            ax_logo = plt.subplot2grid((3,heat_cols), (0,logo_start), colspan=(logo_end-logo_start))
            ax_sad = plt.subplot2grid((3,heat_cols), (1,sad_start), colspan=(sad_end-sad_start))
            ax_heat = plt.subplot2grid((3,heat_cols), (2,0), colspan=heat_cols)

            # print a WebLogo of the sequence
            vlim = max(options.min_limit, abs(minmax_matrix).max())
            if options.gain_height:
                seq_heights = 0.25 + 1.75/vlim*(abs(minmax_matrix).max(axis=0))
            else:
                seq_heights = 0.25 + 1.75/vlim*(-minmax_matrix[0])
            logo_eps = '%s/%s_c%d_seq.eps' % (options.out_dir, header_filename(header), ci)
            seq_logo(seq, seq_heights, logo_eps)

            # add to figure
            logo_png = '%s.png' % logo_eps[:-4]
            logo_cmd = 'convert -density 300 %s %s' % (logo_eps, logo_png)
            if subprocess.call(logo_cmd, shell=True):
                message('Error running convert', 'error')
            logo = Image.open(logo_png)
            ax_logo.imshow(logo)
            ax_logo.set_axis_off()

            # plot loss and gain SAD scores
            ax_sad.plot(-minmax_matrix[0], c=rdbu[0], label='loss', linewidth=1)
            ax_sad.plot(minmax_matrix[1], c=rdbu[-1], label='gain', linewidth=1)
            ax_sad.set_xlim(0,minmax_matrix.shape[1])
            ax_sad.legend()
            # ax_sad.grid(True, linestyle=':')
            for axis in ['top','bottom','left','right']:
                ax_sad.spines[axis].set_linewidth(0.5)

            # plot real-normalized scores
            vlim = max(options.min_limit, abs(norm_matrix).max())
            sns.heatmap(norm_matrix, linewidths=0, cmap='RdBu_r', vmin=-vlim, vmax=vlim, xticklabels=False, ax=ax_heat)
            ax_heat.yaxis.set_ticklabels('TGCA', rotation='horizontal') # , size=10)

            # save final figure
            plt.tight_layout()
            plt.savefig('%s/%s_c%d_heat.pdf' % (options.out_dir,header.replace(':','_'), ci), dpi=300)
            plt.close()


        #################################################################
        # print table of nt variability for each cell
        #################################################################
        for ci in range(seq_mod_preds.shape[3]):
            seq_mod_preds_cell = seq_mod_preds[si,:,:,ci]
            real_pred_cell = get_real_pred(seq_mod_preds_cell, seq)

            min_scores = seq_mod_preds_cell.min(axis=0)
            max_scores = seq_mod_preds_cell.max(axis=0)

            loss_matrix = real_pred_cell - seq_mod_preds_cell.min(axis=0)
            gain_matrix = seq_mod_preds_cell.max(axis=0) - real_pred_cell

            for pos in range(seq_mod_preds_cell.shape[1]):
                cols = [header, delta_start+pos, ci, loss_matrix[pos], gain_matrix[pos]]
                print >> table_out, '\t'.join([str(c) for c in cols])

    table_out.close()


def header_filename(header):
    ''' Revise the FASTA header to be ba better filename '''
    # no colons
    header = header.replace(':','_')

    # no parentheses
    header = header.replace('(','_')
    header = header.replace(')','')

    return header


def get_real_pred(seq_mod_preds, seq):
    ''' Return the real sequence prediction from the modified prediction matrix '''
    si = 0
    while seq[si] == 'N':
        si += 1

    if seq[si] == 'A':
        real_pred = seq_mod_preds[0,si]
    elif seq[si] == 'C':
        real_pred = seq_mod_preds[1,si]
    elif seq[si] == 'G':
        real_pred = seq_mod_preds[2,si]
    else:
        real_pred = seq_mod_preds[3,si]

    return real_pred


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
