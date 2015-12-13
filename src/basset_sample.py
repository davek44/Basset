#!/usr/bin/env python
from optparse import OptionParser
import gzip
import random
import sys

################################################################################
# basset_sample.py
#
# Sample sequences from an existing dataset of sequences as BED file and
# activity table.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <db_bed> <db_act_file> <sample_seqs> <output_prefix>'
    parser = OptionParser(usage)
    parser.add_option('-s', dest='seed', default=1, type='float', help='Random number generator seed [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 4:
    	parser.error('Must provide database BED and activity table and output prefix')
    else:
        bed_file = args[0]
        act_file = args[1]
        sample_seqs = int(args[2])
        out_pre = args[3]

    random.seed(options.seed)

    ############################################################
    # process BED
    ############################################################
    reservoir = ['']*sample_seqs

    bed_in = open(bed_file)

    # initial fill
    i = 0
    while i < sample_seqs:
    	reservoir[i] = bed_in.readline()
    	i += 1

    # sample
    for line in bed_in:
    	j = random.randint(0, i+1)
    	if j < sample_seqs:
    		reservoir[j] = line
    	i += 1

    bed_in.close()

    # print
    bed_out = open('%s.bed' % out_pre, 'w')
    for r in range(len(reservoir)):
        print >> bed_out, reservoir[r],

    bed_out.close()

    ############################################################
    # process activity table
    ############################################################
    # hash sampled headers
    reservoir_headers = set()
    for line in reservoir:
        a = line.rstrip().split('\t')
        chrom = a[0]
        start = a[1]
        end = a[2]
        strand = a[5]
        header = '%s:%s-%s(%s)' % (chrom,start,end,strand)
        reservoir_headers.add(header)

    act_out = open('%s_act.txt' % out_pre, 'w')
    if act_file[-3:] == '.gz':
        act_in = gzip.open(act_file)
    else:
        act_in = open(act_file)

    # print header
    print >> act_out, act_in.readline(),

    # filter activity table
    for line in act_in:
        a = line.split('\t')
        if a[0] in reservoir_headers:
            print >> act_out, line,

    act_in.close()
    act_out.close()



################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
