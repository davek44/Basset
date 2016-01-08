#!/usr/bin/env python
from optparse import OptionParser
import os
import subprocess
import sys
from src.util.utils import message

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
    parser.add_option('--no-cuda', dest='no_cuda', default=False, action='store_true', help='Do not install packages only required for CUDA')
    parser.add_option('-w', dest='warn_on_error', default=False, action='store_true', help='Print a warning, rather than exit, if a dependency cannot be installed [Default: %default]')
    (options,args) = parser.parse_args()

    # confirm luarocks
    message('Checking for luarocks')
    luarocks_which = subprocess.check_output('which luarocks', shell=True)
    if luarocks_which == '':
        message('Please install Torch7 first', 'error')

    ############################################################
    # luarocks database
    ############################################################

    # install luafilesystem
    message('Installing luafilesystem')
    cmd = 'luarocks install luafilesystem'
    if subprocess.call(cmd, shell=True) and not options.warn_on_error:
        sys.exit(1)

    # install dpnn
    message('Installing dpnn')
    cmd = 'luarocks install dpnn'
    if subprocess.call(cmd, shell=True) and not options.warn_on_error:
        sys.exit(1)

    # install inn
    if not options.no_cuda:
        message('Installing inn')
        cmd = 'luarocks install inn'
        if subprocess.call(cmd, shell=True) and not options.warn_on_error:
            sys.exit(1)

    ############################################################
    # luarocks from github
    ############################################################

    os.chdir('src')

    # install torch-hdf5
    cmd = 'git clone https://github.com/deepmind/torch-hdf5.git'
    message('Installing torch-hdf5')
    if subprocess.call(cmd, shell=True) and not options.warn_on_error:
        sys.exit(1)

    os.chdir('torch-hdf5')

    cmd = 'luarocks make'
    message('Installing make')
    if subprocess.call(cmd, shell=True) and not options.warn_on_error:
        sys.exit(1)

    os.chdir('..')

    os.chdir('..')


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
