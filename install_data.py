#!/usr/bin/env python
from optparse import OptionParser
import glob
import os
import subprocess
import sys
from src.util.utils import message

################################################################################
# install_data.py
#
# Download and arrange pre-trained models and data.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    parser.add_option('-r', dest='restart', default=False, action='store_true', help='Do not overwrite existing files, as if restarting an aborted installation [Default: %default]')
    parser.add_option('-w', dest='warn_on_error', default=False, action='store_true', help='Print a warning, rather than exit, if a dependency cannot be installed [Default: %default]')
    (options,args) = parser.parse_args()

    os.chdir('data')

    ############################################################
    # download pre-trained model
    ############################################################
    os.chdir('models')

    if not options.restart or not os.path.isfile('pretrained_model.th'):
        message('Downloading pre-trained model.')

        cmd = 'wget https://www.dropbox.com/s/rguytuztemctkf8/pretrained_model.th.gz'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error:
            sys.exit(1)

        cmd = 'gunzip pretrained_model.th.gz'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error:
            sys.exit(1)

    os.chdir('..')


    ############################################################
    # download human genome
    ############################################################
    os.chdir('genomes')

    if not options.restart or not os.path.isfile('hg19.fa'):
        message('Downloading hg19 FASTA from UCSC. If you already have it, CTL-C to place a sym link in the genomes directory named hg19.fa')

        # download hg19
        cmd = 'wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz -O chromFa.tar.gz'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error:
            sys.exit(1)

        # un-tar
        cmd = 'tar -xzvf chromFa.tar.gz'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error:
            sys.exit(1)

        # cat
        cmd = 'cat chr?.fa chr??.fa > hg19.fa'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error:
            sys.exit(1)

        # clean up
        os.remove('chromFa.tar.gz')
        for chrom_fa in glob.glob('chr*.fa'):
            os.remove(chrom_fa)

    if not options.restart or not os.path.isfile('hg19.fa.fai'):
        cmd = 'samtools faidx hg19.fa'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error:
            sys.exit(1)

    os.chdir('..')


    ############################################################
    # download and prepare public data
    ############################################################
    if not options.restart or not os.path.isfile('encode_roadmap.h5'):
        message('Downloading and preparing public data')
        cmd = 'wget https://www.dropbox.com/s/h1cqokbr8vjj5wc/encode_roadmap.bed.gz'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error: sys.exit(1)
        cmd = 'gunzip encode_roadmap.bed.gz'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error: sys.exit(1)

        cmd = 'wget https://www.dropbox.com/s/8g3kc0ai9ir5d15/encode_roadmap_act.txt.gz'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error: sys.exit(1)
        cmd = 'gunzip encode_roadmap_act.txt.gz'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error: sys.exit(1)

        '''
        # download and arrange available data
        cmd = './get_dnase.sh'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error: sys.exit(1)

        # preprocess
        cmd = 'preprocess_features.py -y -m 200 -s 600 -o encode_roadmap -c human.hg19.genome sample_beds.txt'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error: sys.exit(1)
        '''

        # make a FASTA file
        cmd = 'bedtools getfasta -fi genomes/hg19.fa -bed encode_roadmap.bed -s -fo encode_roadmap.fa'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error: sys.exit(1)

        # make an HDF5 file
        cmd = 'seq_hdf5.py -c -r -t 71886 -v 70000 encode_roadmap.fa encode_roadmap_act.txt encode_roadmap.h5'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error: sys.exit(1)



################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
