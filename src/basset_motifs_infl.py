#!/usr/bin/env python
from optparse import OptionParser
import os
import random
import sys

import numpy as np
import pandas as pd
from sklearn import preprocessing
import h5py
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
from util.utils import message

################################################################################
# basset_motifs_infl.py
#
# Collect statistics and make plots to explore the influence of filters in the
# first convolution layer of the given model using the given sequences.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <model_file> <test_hdf5_file>'
    parser = OptionParser(usage)
    parser.add_option('-b', dest='batch_size', default=1000, type='int', help='Batch size (affects memory usage) [Default: %default]')
    parser.add_option('-d', dest='model_hdf5_file', default=None, help='Pre-computed model output as HDF5.')
    parser.add_option('-i', dest='informative_only', default=False, action='store_true', help='Plot informative filters only [Default: %default]')
    parser.add_option('-m', dest='motifs_file')
    parser.add_option('-n', dest='norm_targets', default=False, action='store_true', help='Use the norm of the target influences as the primary influence measure [Default: %default]')
    parser.add_option('-o', dest='out_dir', default='.')
    parser.add_option('--subset', dest='subset_file', default=None, help='Subset targets to the ones in this file')
    parser.add_option('-s', dest='sample', default=None, type='int', help='Sample sequences from the test set [Default:%default]')
    parser.add_option('--seqs', dest='seqs', default=False, action='store_true', help='Output sequence-specific influence [Default: %default]')
    parser.add_option('-t', dest='targets_file', default=None, help='File specifying target indexes and labels in table format')
    parser.add_option('--width', dest='heat_width', default=10, type='float')
    parser.add_option('--height', dest='heat_height', default=20, type='float')
    parser.add_option('--font', dest='heat_font', default=0.4, type='float')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide Basset model file and test data in HDF5 format.')
    else:
        model_file = args[0]
        test_hdf5_file = args[1]

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    #################################################################
    # load data
    #################################################################
    # load sequences
    test_hdf5_in = h5py.File(test_hdf5_file, 'r')
    seq_vecs = np.array(test_hdf5_in['test_in'])
    seq_targets = np.array(test_hdf5_in['test_out'])
    seq_headers = np.array(test_hdf5_in['test_headers'])
    test_hdf5_in.close()

    # name the targets
    target_names = name_targets(seq_targets.shape[1], options.targets_file)

    if options.subset_file:
        target_subset = set([line.rstrip() for line in open(options.subset_file)])

    # get additional motif information
    df_motifs = None
    if options.motifs_file:
        df_motifs = pd.read_table(options.motifs_file, delim_whitespace=True)


    #################################################################
    # sample
    #################################################################
    if options.sample is not None:
        # choose sampled indexes
        sample_i = np.array(random.sample(xrange(seq_vecs.shape[0]), options.sample))

        # filter
        seq_vecs = seq_vecs[sample_i]
        seq_targets = seq_targets[sample_i]
        seq_headers = seq_headers[sample_i]

        # create a new HDF5 file
        sample_hdf5_file = '%s/sample.h5' % options.out_dir
        sample_hdf5_out = h5py.File(sample_hdf5_file, 'w')
        sample_hdf5_out.create_dataset('test_in', data=seq_vecs)
        sample_hdf5_out.create_dataset('test_out', data=seq_targets)
        sample_hdf5_out.close()

        # update test HDF5
        test_hdf5_file = sample_hdf5_file


    #################################################################
    # Torch predict
    #################################################################
    if options.model_hdf5_file is None:
        torch_opts = ''
        if options.seqs:
            torch_opts += '-seqs'

        options.model_hdf5_file = '%s/model_out.h5' % options.out_dir
        torch_cmd = 'basset_motifs_infl.lua -batch_size %d %s %s %s %s' % (options.batch_size, torch_opts, model_file, test_hdf5_file, options.model_hdf5_file)
        if subprocess.call(torch_cmd, shell=True):
            message('Error running basset_motifs_infl.lua', 'error')

    # load model output
    model_hdf5_in = h5py.File(options.model_hdf5_file, 'r')
    filter_means = np.array(model_hdf5_in['filter_means'])
    filter_stds = np.array(model_hdf5_in['filter_stds'])
    filter_infl = np.array(model_hdf5_in['filter_infl'])
    filter_infl_targets = np.array(model_hdf5_in['filter_infl_targets'])
    if options.seqs:
        seq_filter_targets = np.array(model_hdf5_in['seq_filter_targets'])
    model_hdf5_in.close()


    #############################################################
    # use target-based influence
    #############################################################
    if options.norm_targets:
        # save the loss-based influence
        filter_infl_loss = np.array(filter_infl, copy=True)

        # set to the target-based influence
        for fi in range(filter_infl_targets.shape[0]):
            filter_infl[fi] = np.mean(filter_infl_targets[fi]**2)

        # print to a table
        tnorm_out = open('%s/loss_target.txt' % options.out_dir, 'w')
        for fi in range(len(filter_infl)):
            cols = (fi, filter_infl_loss[fi], filter_infl[fi])
            print >> tnorm_out, '%3d  %7.4f  %7.4f' % cols
        tnorm_out.close()

        # compare the two
        xmin, xmax = coord_range(filter_infl_loss, buf_pct=0.1)
        ymin, ymax = coord_range(filter_infl, buf_pct=0.1)

        sns.set(style='ticks', font_scale=1)
        plt.figure()
        g = sns.jointplot(x=filter_infl_loss, y=filter_infl, color='black', joint_kws={'alpha':0.7})
        ax = g.ax_joint
        ax.set_xlim(xmin, xmax)
        ax.set_xlabel('loss-based influence')
        ax.xaxis.label.set_fontsize(18)
        map(lambda xl: xl.set_fontsize(15), ax.get_xticklabels())
        ax.set_ylim(ymin, ymax)
        ax.set_ylabel('target-based influence')
        ax.yaxis.label.set_fontsize(18)
        map(lambda yl: yl.set_fontsize(15), ax.get_yticklabels())
        ax.grid(True, linestyle=':')
        plt.tight_layout(w_pad=0, h_pad=0)
        plt.savefig('%s/loss_target.pdf' % options.out_dir)
        plt.close()


    #############################################################
    # print filter influence table
    #############################################################
    table_out = open('%s/table.txt' % options.out_dir, 'w')
    for i in range(len(filter_infl)):
        if df_motifs is not None:
            cols = (i, filter_infl[i], filter_means[i], filter_stds[i], df_motifs.ic.iloc[i], df_motifs.annotation.iloc[i])
            print >> table_out, '%3d  %7.4f  %6.4f  %6.3f  %4.1f  %8s' % cols
        else:
            cols = (i, filter_infl[i], filter_means[i], filter_stds[i])
            print >> table_out, '%3d  %7.4f  %6.4f  %6.3f' % cols
    table_out.close()


    #################################################################
    # plot filter influence
    #################################################################
    sb_blue = sns.color_palette('deep')[0]
    sns.set(style='ticks', font_scale=1)
    ymin, ymax = coord_range(filter_infl, buf_pct=0.1)

    if options.motifs_file:
        nonzero = np.array(df_motifs.ic > 0)
        xmin, xmax = coord_range(df_motifs.ic.loc[nonzero])
        plt.figure()
        g = sns.jointplot(x=np.array(df_motifs.ic.loc[nonzero]), y=filter_infl[nonzero], color='black', stat_func=None, joint_kws={'alpha':0.8})
        ax = g.ax_joint
        ax.set_xlim(xmin, xmax)
        ax.set_xlabel('Information content')
        ax.xaxis.label.set_fontsize(18)
        map(lambda xl: xl.set_fontsize(15), ax.get_xticklabels())
        ax.set_ylim(ymin, ymax)
        ax.set_ylabel('Influence')
        ax.yaxis.label.set_fontsize(18)
        map(lambda yl: yl.set_fontsize(15), ax.get_yticklabels())
        # ax.grid(True, linestyle=':')
        plt.tight_layout(w_pad=0, h_pad=0)
        plt.savefig('%s/ic_infl.pdf' % options.out_dir)
        plt.close()


    #############################################################
    # prep for cell-specific analyses
    #############################################################
    filter_names = name_filters(len(filter_infl), df_motifs)

    # construct a panda data frame of the target influences
    df_ft = pd.DataFrame(filter_infl_targets, index=filter_names, columns=target_names)

    # print filter influence per target table
    table_out = open('%s/table_target.txt' % options.out_dir, 'w')
    for i in range(df_ft.shape[0]):
        for ti in range(len(target_names)):
            cols = (i, ti, target_names[ti], df_ft.iloc[i,ti])
            print >> table_out, '%-3d  %3d  %20s  %7.4f' % cols
    table_out.close()

    # print sequence-specific filter influence per target table
    if options.seqs:
        table_out = open('%s/table_seqs.txt' % options.out_dir, 'w')
        for si in range(seq_filter_targets.shape[0]):
            for fi in range(seq_filter_targets.shape[1]):
                for ti in range(seq_filter_targets.shape[2]):
                    cols = (seq_headers[si], fi, ti, seq_filter_targets[si][fi][ti])
                    print >> table_out, '%-25s  %3d  %3d  %7.4f' % cols
        table_out.close()

    # use only high information filters
    if options.informative_only and df_motifs is not None:
        df_ft = df_ft[df_moitfs.ic > 6]
    elif df_ft.shape[1] >= 10:
        df_ft_stds = df_ft.std(axis=1)
        df_ft = df_ft[df_ft_stds > 0]

    #############################################################
    # plot filter influence per cell heatmaps
    #############################################################
    # subset targets before plotting
    if options.subset_file:
        subset_mask = df_ft.columns.isin(target_subset)
        df_ft_sub = df_ft.loc[:,subset_mask]

        plot_infl_heatmaps(df_ft_sub, options.out_dir, options.heat_width, options.heat_height, options.heat_font)

    # plot all cells
    else:
        plot_infl_heatmaps(df_ft, options.out_dir, options.heat_width, options.heat_height, options.heat_font)


def coord_range(nums, buf_pct=0.05):
    ''' Determine a nice buffered axis range from a list/array of numbers '''
    nmin = min(nums)
    nmax = max(nums)
    spread = nmax-nmin
    buf = buf_pct*spread
    return nmin-buf, nmax+buf


def name_filters(num_filters, df_motifs):
    ''' Name the filters using Tomtom matches.

    Attrs:
        num_filters (int) : total number of filters
        df_motifs (DataFrame) : DataFrame with filter annotations

    Returns:
        filter_names [str] :
    '''
    # name by number
    filter_names = ['f%d'%fi for fi in range(num_filters)]

    # name by protein
    if df_motifs is not None:
        for fi in range(num_filters):
            ann = df_motifs.annotation.iloc[fi]
            if ann != '.':
                filter_names[fi] += '_%s' % ann

    return np.array(filter_names)


def name_targets(num_targets, targets_file):
    ''' Name the targets using a file of names.

    Attrs:
        num_targets (int)
        targets_file (str)
    Returns:
        target_names [str]
    '''
    if targets_file == None:
        target_names = ['t%d' % ti for ti in range(num_targets)]
    else:
        target_names = [line.split()[1] for line in open(targets_file)]
    return target_names


def plot_heat(mat, out_pdf, width, height):
    ''' Plot a single filter influence heat map'''

    if mat.shape[1] > 2:
        plt.figure()
        g = sns.clustermap(mat, figsize=(width,height))

        for tick in g.ax_heatmap.get_xticklabels():
            tick.set_rotation(-45)
            tick.set_horizontalalignment('left')
    else:
        plt.figure(figsize=(width,height))
        g = sns.heatmap(mat)

    plt.savefig(out_pdf)
    plt.close()


def plot_infl_heatmaps(unit_target_deltas, out_dir, width=10, height=20, font=0.4):
    ''' Plot a variety of heatmaps about the filter influence per cell data frame '''

    sns.set(style='white', font_scale=font)

    # plot influences per cell
    plot_heat(unit_target_deltas, '%s/infl_target.pdf' % out_dir, width, height)

    # normalize per cell
    utd_z = preprocessing.scale(unit_target_deltas, axis=1)
    unit_target_deltas_norm = pd.DataFrame(utd_z, index=unit_target_deltas.index, columns=unit_target_deltas.columns)

    plot_heat(unit_target_deltas_norm, '%s/infl_target_norm.pdf' % out_dir, width, height)

    # use only annotated filters
    annotated = np.array([label.find('_')!=-1 for label in unit_target_deltas_norm.index])
    unit_target_deltas_ann = unit_target_deltas_norm[annotated]

    if unit_target_deltas_ann.shape[0] > 0:
        plot_heat(unit_target_deltas_ann, '%s/infl_target_ann.pdf' % out_dir, width, height)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
