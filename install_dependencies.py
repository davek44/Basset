#!/usr/bin/env python  
from optparse import OptionParser
import os
import subprocess
import sys

################################################################################
# install_dependencies.py  
#
# Download and install Basset dependencies.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()

    # confirm luarocks
    luarocks_which = subprocess.check_output('which luarocks', shell=True)
    if luarocks_which == '':
        print >> sys.stderr, 'Please install Torch7 first.'
        exit(1)

    ############################################################
    # luarocks database
    ############################################################

    # install luafilesystem
    cmd = 'luarocks install luafilesystem'
    subprocess.call(cmd, shell=True)

    # install dpnn
    cmd = 'luarocks install dpnn'
    subprocess.call(cmd, shell=True)

    # install inn
    cmd = 'luarocks install inn'
    subprocess.call(cmd, shell=True)

    # install dp
    cmd = 'luarocks install dp'
    subprocess.call(cmd, shell=True)


    ############################################################
    # luarocks from github
    ############################################################

    os.chdir('src')

    # install torch-hdf5
    cmd = 'git clone https://github.com/davek44/torch-hdf5.git'
    subprocess.call(cmd, shell=True)

    os.chdir('torch-hdf5')

    cmd = 'luarocks make'
    subprocess.call(cmd, shell=True)

    os.chdir('..')

    os.chdir('..')


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
