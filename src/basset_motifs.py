#!/usr/bin/env python
from optparse import OptionParser
import copy, os, pdb, random, shutil, time

import h5py
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import seaborn as sns
from sklearn import preprocessing
import subprocess
import dna_io
from util.utils import message

################################################################################
# basset_motifs.py
#
# Collect statistics and make plots to explore the first convolution layer
# of the given model using the given sequences.
################################################################################

weblogo_opts = '-X NO -Y NO --errorbars NO --fineprint ""'
weblogo_opts += ' -C "#CB2026" A A'
weblogo_opts += ' -C "#34459C" C C'
weblogo_opts += ' -C "#FBB116" G G'
weblogo_opts += ' -C "#0C8040" T T'

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <model_file> <test_hdf5_file>'
    parser = OptionParser(usage)
    parser.add_option('-d', dest='model_hdf5_file', default=None, help='Pre-computed model output as HDF5.')
    parser.add_option('-o', dest='out_dir', default='.')
    parser.add_option('-m', dest='meme_db', default='%s/data/motifs/Homo_sapiens.meme' % os.environ['BASSETDIR'], help='MEME database used to annotate motifs')
    parser.add_option('-s', dest='sample', default=None, type='int', help='Sample sequences from the test set [Default:%default]')
    parser.add_option('-t', dest='trim_filters', default=False, action='store_true', help='Trim uninformative positions off the filter ends [Default: %default]')
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
    try:
        target_names = list(test_hdf5_in['target_labels'])
    except KeyError:
        target_names = ['t%d'%ti for ti in range(seq_targets.shape[1])]
    test_hdf5_in.close()


    #################################################################
    # sample
    #################################################################
    if options.sample is not None:
        # choose sampled indexes
        sample_i = np.array(random.sample(xrange(seq_vecs.shape[0]), options.sample))

        # filter
        seq_vecs = seq_vecs[sample_i]
        seq_targets = seq_targets[sample_i]

        # create a new HDF5 file
        sample_hdf5_file = '%s/sample.h5' % options.out_dir
        sample_hdf5_out = h5py.File(sample_hdf5_file, 'w')
        sample_hdf5_out.create_dataset('test_in', data=seq_vecs)
        sample_hdf5_out.close()

        # update test HDF5
        test_hdf5_file = sample_hdf5_file

    # convert to letters
    seqs = dna_io.vecs2dna(seq_vecs)


    #################################################################
    # Torch predict
    #################################################################
    if options.model_hdf5_file is None:
        options.model_hdf5_file = '%s/model_out.h5' % options.out_dir
        torch_cmd = 'basset_motifs_predict.lua %s %s %s' % (model_file, test_hdf5_file, options.model_hdf5_file)
        if subprocess.call(torch_cmd, shell=True):
            message('Error running basset_motifs_predict.lua', 'error')

    # load model output
    model_hdf5_in = h5py.File(options.model_hdf5_file, 'r')
    filter_weights = np.array(model_hdf5_in['weights'])
    filter_outs = np.array(model_hdf5_in['outs'])
    model_hdf5_in.close()

    # store useful variables
    num_filters = filter_weights.shape[0]
    filter_size = filter_weights.shape[2]


    #################################################################
    # individual filter plots
    #################################################################
    # also save information contents
    filters_ic = []
    meme_out = meme_intro('%s/filters_meme.txt'%options.out_dir, seqs)

    for f in range(num_filters):
        print 'Filter %d' % f

        # plot filter parameters as a heatmap
        plot_filter_heat(filter_weights[f,:,:], '%s/filter%d_heat.pdf' % (options.out_dir,f))

        # write possum motif file
        filter_possum(filter_weights[f,:,:], 'filter%d'%f, '%s/filter%d_possum.txt'%(options.out_dir,f), options.trim_filters)

        # plot weblogo of high scoring outputs
        plot_filter_logo(filter_outs[:,f,:], filter_size, seqs, '%s/filter%d_logo'%(options.out_dir,f), maxpct_t=0.5)

        # make a PWM for the filter
        filter_pwm, nsites = make_filter_pwm('%s/filter%d_logo.fa'%(options.out_dir,f))

        if nsites < 10:
            # no information
            filters_ic.append(0)
        else:
            # compute and save information content
            filters_ic.append(info_content(filter_pwm))

            # add to the meme motif file
            meme_add(meme_out, f, filter_pwm, nsites, options.trim_filters)

    meme_out.close()


    #################################################################
    # annotate filters
    #################################################################
    # run tomtom
    tomtom_cmd = 'tomtom -dist pearson -thresh 0.1 -oc %s/tomtom %s/filters_meme.txt %s' % (options.out_dir, options.out_dir, options.meme_db)
    if subprocess.call(tomtom_cmd, shell=True):
        message('Error running tomtom', 'error')

    # read in annotations
    filter_names = name_filters(num_filters, '%s/tomtom/tomtom.txt'%options.out_dir, options.meme_db)


    #################################################################
    # print a table of information
    #################################################################
    table_out = open('%s/table.txt'%options.out_dir, 'w')

    # print header for later panda reading
    header_cols = ('', 'consensus', 'annotation', 'ic', 'mean', 'std')
    print >> table_out, '%3s  %19s  %10s  %5s  %6s  %6s' % header_cols

    for f in range(num_filters):
        # collapse to a consensus motif
        consensus = filter_motif(filter_weights[f,:,:])

        # grab annotation
        annotation = '.'
        name_pieces = filter_names[f].split('_')
        if len(name_pieces) > 1:
            annotation = name_pieces[1]

        # plot density of filter output scores
        fmean, fstd = plot_score_density(np.ravel(filter_outs[:,f,:]), '%s/filter%d_dens.pdf' % (options.out_dir,f))

        row_cols = (f, consensus, annotation, filters_ic[f], fmean, fstd)
        print >> table_out, '%-3d  %19s  %10s  %5.2f  %6.4f  %6.4f' % row_cols

    table_out.close()


    #################################################################
    # global filter plots
    #################################################################
    # plot filter-sequence heatmap
    plot_filter_seq_heat(filter_outs, '%s/filter_seqs.pdf'%options.out_dir)

    # plot filter-segment heatmap
    plot_filter_seg_heat(filter_outs, '%s/filter_segs.pdf'%options.out_dir)
    plot_filter_seg_heat(filter_outs, '%s/filter_segs_raw.pdf'%options.out_dir, whiten=False)

    # plot filter-target correlation heatmap
    plot_target_corr(filter_outs, seq_targets, filter_names, target_names, '%s/filter_target_cors_mean.pdf'%options.out_dir, 'mean')
    plot_target_corr(filter_outs, seq_targets, filter_names, target_names, '%s/filter_target_cors_max.pdf'%options.out_dir, 'max')


def get_motif_proteins(meme_db_file):
    ''' Hash motif_id's to protein names using the MEME DB file '''
    motif_protein = {}
    for line in open(meme_db_file):
        a = line.split()
        if len(a) > 0 and a[0] == 'MOTIF':
            if a[2][0] == '(':
                motif_protein[a[1]] = a[2][1:a[2].find(')')]
            else:
                motif_protein[a[1]] = a[2]
    return motif_protein


def info_content(pwm, transpose=False):
    ''' Compute PWM information content '''
    pseudoc = 1e-9

    if transpose:
        pwm = np.transpose(pwm)

    ic = 0
    for i in range(pwm.shape[0]):
        for j in range(4):
            ic += 0.5 + pwm[i][j]*np.log2(pseudoc+pwm[i][j])

    return ic


def make_filter_pwm(filter_fasta):
    ''' Make a PWM for this filter from its top hits '''

    nts = {'A':0, 'C':1, 'G':2, 'T':3}
    pwm_counts = []
    nsites = 4 # pseudocounts
    for line in open(filter_fasta):
        if line[0] != '>':
            seq = line.rstrip()
            nsites += 1
            if len(pwm_counts) == 0:
                # initialize with the length
                for i in range(len(seq)):
                    pwm_counts.append(np.array([1.0]*4))

            # count
            for i in range(len(seq)):
                try:
                    pwm_counts[i][nts[seq[i]]] += 1
                except KeyError:
                    pwm_counts[i] += np.array([0.25]*4)

    # normalize
    pwm_freqs = []
    for i in range(len(pwm_counts)):
        pwm_freqs.append([pwm_counts[i][j]/float(nsites) for j in range(4)])

    return np.array(pwm_freqs), nsites-4


def meme_add(meme_out, f, filter_pwm, nsites, trim_filters=False):
    ''' Print a filter to the growing MEME file

    Attrs:
        meme_out : open file
        f (int) : filter index #
        filter_pwm (array) : filter PWM array
        nsites (int) : number of filter sites
    '''
    if not trim_filters:
        ic_start = 0
        ic_end = filter_pwm.shape[0]-1
    else:
        ic_t = 0.2

        # trim PWM of uninformative prefix
        ic_start = 0
        while ic_start < filter_pwm.shape[0] and info_content(filter_pwm[ic_start:ic_start+1]) < ic_t:
            ic_start += 1

        # trim PWM of uninformative suffix
        ic_end = filter_pwm.shape[0]-1
        while ic_end >= 0 and info_content(filter_pwm[ic_end:ic_end+1]) < ic_t:
            ic_end -= 1

    if ic_start < ic_end:
        print >> meme_out, 'MOTIF filter%d' % f
        print >> meme_out, 'letter-probability matrix: alength= 4 w= %d nsites= %d' % (ic_end-ic_start+1, nsites)

        for i in range(ic_start, ic_end+1):
            print >> meme_out, '%.4f %.4f %.4f %.4f' % tuple(filter_pwm[i])
        print >> meme_out, ''


def meme_intro(meme_file, seqs):
    ''' Open MEME motif format file and print intro

    Attrs:
        meme_file (str) : filename
        seqs [str] : list of strings for obtaining background freqs

    Returns:
        mem_out : open MEME file
    '''
    nts = {'A':0, 'C':1, 'G':2, 'T':3}

    # count
    nt_counts = [1]*4
    for i in range(len(seqs)):
        for nt in seqs[i]:
            try:
                nt_counts[nts[nt]] += 1
            except KeyError:
                pass

    # normalize
    nt_sum = float(sum(nt_counts))
    nt_freqs = [nt_counts[i]/nt_sum for i in range(4)]

    # open file for writing
    meme_out = open(meme_file, 'w')

    # print intro material
    print >> meme_out, 'MEME version 4'
    print >> meme_out, ''
    print >> meme_out, 'ALPHABET= ACGT'
    print >> meme_out, ''
    print >> meme_out, 'Background letter frequencies:'
    print >> meme_out, 'A %.4f C %.4f G %.4f T %.4f' % tuple(nt_freqs)
    print >> meme_out, ''

    return meme_out


def name_filters(num_filters, tomtom_file, meme_db_file):
    ''' Name the filters using Tomtom matches.

    Attrs:
        num_filters (int) : total number of filters
        tomtom_file (str) : filename of Tomtom output table.
        meme_db_file (str) : filename of MEME db

    Returns:
        filter_names [str] :
    '''
    # name by number
    filter_names = ['f%d'%fi for fi in range(num_filters)]

    # name by protein
    if tomtom_file is not None and meme_db_file is not None:
        motif_protein = get_motif_proteins(meme_db_file)

        # hash motifs and q-value's by filter
        filter_motifs = {}

        tt_in = open(tomtom_file)
        tt_in.readline()
        for line in tt_in:
            a = line.split()
            fi = int(a[0][6:])
            motif_id = a[1]
            qval = float(a[5])

            filter_motifs.setdefault(fi,[]).append((qval,motif_id))

        tt_in.close()

        # assign filter's best match
        for fi in filter_motifs:
            top_motif = sorted(filter_motifs[fi])[0][1]
            filter_names[fi] += '_%s' % motif_protein[top_motif]

    return np.array(filter_names)


################################################################################
# plot_target_corr
#
# Plot a clustered heatmap of correlations between filter activations and
# targets.
#
# Input
#  filter_outs:
#  filter_names:
#  target_names:
#  out_pdf:
################################################################################
def plot_target_corr(filter_outs, seq_targets, filter_names, target_names, out_pdf, seq_op='mean'):
    num_seqs = filter_outs.shape[0]
    num_targets = len(target_names)

    if seq_op == 'mean':
        filter_outs_seq = filter_outs.mean(axis=2)
    else:
        filter_outs_seq = filter_outs.max(axis=2)

    # std is sequence by filter.
    filter_seqs_std = filter_outs_seq.std(axis=0)
    filter_outs_seq = filter_outs_seq[:,filter_seqs_std > 0]
    filter_names_live = filter_names[filter_seqs_std > 0]

    filter_target_cors = np.zeros((len(filter_names_live),num_targets))
    for fi in range(len(filter_names_live)):
        for ti in range(num_targets):
            cor, p = spearmanr(filter_outs_seq[:,fi], seq_targets[:num_seqs,ti])
            filter_target_cors[fi,ti] = cor

    cor_df = pd.DataFrame(filter_target_cors, index=filter_names_live, columns=target_names)

    sns.set(font_scale=0.3)
    plt.figure()
    sns.clustermap(cor_df, cmap='BrBG', center=0, figsize=(8,10))
    plt.savefig(out_pdf)
    plt.close()


################################################################################
# plot_filter_seq_heat
#
# Plot a clustered heatmap of filter activations in
#
# Input
#  param_matrix: np.array of the filter's parameter matrix
#  out_pdf:
################################################################################
def plot_filter_seq_heat(filter_outs, out_pdf, whiten=True, drop_dead=True):
    # compute filter output means per sequence
    filter_seqs = filter_outs.mean(axis=2)

    # whiten
    if whiten:
        filter_seqs = preprocessing.scale(filter_seqs)

    # transpose
    filter_seqs = np.transpose(filter_seqs)

    if drop_dead:
        filter_stds = filter_seqs.std(axis=1)
        filter_seqs = filter_seqs[filter_stds > 0]

    # downsample sequences
    seqs_i = np.random.randint(0, filter_seqs.shape[1], 500)

    hmin = np.percentile(filter_seqs[:,seqs_i], 0.1)
    hmax = np.percentile(filter_seqs[:,seqs_i], 99.9)

    sns.set(font_scale=0.3)

    plt.figure()
    sns.clustermap(filter_seqs[:,seqs_i], row_cluster=True, col_cluster=True, linewidths=0, xticklabels=False, vmin=hmin, vmax=hmax)
    plt.savefig(out_pdf)
    #out_png = out_pdf[:-2] + 'ng'
    #plt.savefig(out_png, dpi=300)
    plt.close()


################################################################################
# plot_filter_seq_heat
#
# Plot a clustered heatmap of filter activations in sequence segments.
#
# Mean doesn't work well for the smaller segments for some reason, but taking
# the max looks OK. Still, similar motifs don't cluster quite as well as you
# might expect.
#
# Input
#  filter_outs
################################################################################
def plot_filter_seg_heat(filter_outs, out_pdf, whiten=True, drop_dead=True):
    b = filter_outs.shape[0]
    f = filter_outs.shape[1]
    l = filter_outs.shape[2]

    s = 5
    while l/float(s) - (l/s) > 0:
        s += 1
    print '%d segments of length %d' % (s,l/s)

    # split into multiple segments
    filter_outs_seg = np.reshape(filter_outs, (b, f, s, l/s))

    # mean across the segments
    filter_outs_mean = filter_outs_seg.max(axis=3)

    # break each segment into a new instance
    filter_seqs = np.reshape(np.swapaxes(filter_outs_mean, 2, 1), (s*b, f))

    # whiten
    if whiten:
        filter_seqs = preprocessing.scale(filter_seqs)

    # transpose
    filter_seqs = np.transpose(filter_seqs)

    if drop_dead:
        filter_stds = filter_seqs.std(axis=1)
        filter_seqs = filter_seqs[filter_stds > 0]

    # downsample sequences
    seqs_i = np.random.randint(0, filter_seqs.shape[1], 500)

    hmin = np.percentile(filter_seqs[:,seqs_i], 0.1)
    hmax = np.percentile(filter_seqs[:,seqs_i], 99.9)

    sns.set(font_scale=0.3)
    if whiten:
        dist = 'euclidean'
    else:
        dist = 'cosine'

    plt.figure()
    sns.clustermap(filter_seqs[:,seqs_i], metric=dist, row_cluster=True, col_cluster=True, linewidths=0, xticklabels=False, vmin=hmin, vmax=hmax)
    plt.savefig(out_pdf)
    #out_png = out_pdf[:-2] + 'ng'
    #plt.savefig(out_png, dpi=300)
    plt.close()


################################################################################
# filter_motif
#
# Collapse the filter parameter matrix to a single DNA motif.
#
# Input
#  param_matrix: np.array of the filter's parameter matrix
#  out_pdf:
################################################################################
def filter_motif(param_matrix):
    nts = 'ACGT'

    motif_list = []
    for v in range(param_matrix.shape[1]):
        max_n = 0
        for n in range(1,4):
            if param_matrix[n,v] > param_matrix[max_n,v]:
                max_n = n

        if param_matrix[max_n,v] > 0:
            motif_list.append(nts[max_n])
        else:
            motif_list.append('N')

    return ''.join(motif_list)


################################################################################
# filter_possum
#
# Write a Possum-style motif
#
# Input
#  param_matrix: np.array of the filter's parameter matrix
#  out_pdf:
################################################################################
def filter_possum(param_matrix, motif_id, possum_file, trim_filters=False, mult=200):
    # possible trim
    trim_start = 0
    trim_end = param_matrix.shape[1]-1
    trim_t = 0.3
    if trim_filters:
        # trim PWM of uninformative prefix
        while np.max(param_matrix[:,trim_start]) - np.min(param_matrix[:,trim_start]) < trim_t:
            trim_start += 1

        # trim PWM of uninformative suffix
        while np.max(param_matrix[:,trim_end]) - np.min(param_matrix[:,trim_end]) < trim_t:
            trim_end -= 1

    possum_out = open(possum_file, 'w')
    print >> possum_out, 'BEGIN GROUP'
    print >> possum_out, 'BEGIN FLOAT'
    print >> possum_out, 'ID %s' % motif_id
    print >> possum_out, 'AP DNA'
    print >> possum_out, 'LE %d' % (trim_end+1-trim_start)
    for ci in range(trim_start,trim_end+1):
        print >> possum_out, 'MA %s' % ' '.join(['%.2f'%(mult*n) for n in param_matrix[:,ci]])
    print >> possum_out, 'END'
    print >> possum_out, 'END'

    possum_out.close()


################################################################################
# plot_filter_heat
#
# Plot a heatmap of the filter's parameters.
#
# Input
#  param_matrix: np.array of the filter's parameter matrix
#  out_pdf:
################################################################################
def plot_filter_heat(param_matrix, out_pdf):
    param_range = abs(param_matrix).max()

    plt.figure(figsize=(param_matrix.shape[1], 4))
    sns.heatmap(param_matrix, cmap='PRGn', linewidths=0.2, vmin=-param_range, vmax=param_range)
    ax = plt.gca()
    ax.set_xticklabels(range(1,param_matrix.shape[1]+1))
    ax.set_yticklabels('TGCA', rotation='horizontal', size=10)
    plt.savefig(out_pdf)
    plt.close()


################################################################################
# plot_filter_logo
#
# Plot a weblogo of the filter's occurrences
#
# Input
#  param_matrix: np.array of the filter's parameter matrix
#  out_pdf:
################################################################################
def plot_filter_logo(filter_outs, filter_size, seqs, out_prefix, raw_t=0, maxpct_t=None):
    if maxpct_t:
        all_outs = np.ravel(filter_outs)
        all_outs_mean = all_outs.mean()
        all_outs_norm = all_outs - all_outs_mean
        raw_t = maxpct_t * all_outs_norm.max() + all_outs_mean

    # print fasta file of positive outputs
    filter_fasta_out = open('%s.fa' % out_prefix, 'w')
    filter_count = 0
    for i in range(filter_outs.shape[0]):
        for j in range(filter_outs.shape[1]):
            if filter_outs[i,j] > raw_t:
                kmer = seqs[i][j:j+filter_size]
                print >> filter_fasta_out, '>%d_%d' % (i,j)
                print >> filter_fasta_out, kmer
                filter_count += 1
    filter_fasta_out.close()

    # make weblogo
    if filter_count > 0:
        weblogo_cmd = 'weblogo %s < %s.fa > %s.eps' % (weblogo_opts, out_prefix, out_prefix)
        if subprocess.call(weblogo_cmd, shell=True):
            message('Error running weblogo', 'error')


################################################################################
# plot_score_density
#
# Plot the score density and print to the stats table.
#
# Input
#  param_matrix: np.array of the filter's parameter matrix
#  out_pdf:
################################################################################
def plot_score_density(f_scores, out_pdf):
    sns.set(font_scale=1.3)
    plt.figure()
    sns.distplot(f_scores, kde=False)
    plt.xlabel('ReLU output')
    plt.savefig(out_pdf)
    plt.close()

    return f_scores.mean(), f_scores.std()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
