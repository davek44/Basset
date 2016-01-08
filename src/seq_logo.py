#!/usr/bin/env python
from optparse import OptionParser
import os, re, tempfile
import subprocess
from util.utils import message

################################################################################
# seq_logo
#
# Construct an arbitrary sequence logo using weblogo
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    seq = 'ACGTACGT'
    heights = [1, 1, 2, 2, 1, 1, 1, 1]
    out_eps = 'test_logo.eps'
    seq_logo(seq, heights, out_eps)


def seq_logo(seq, heights, out_eps, weblogo_args=''):
	# print the sequence to a temp fasta file
	fasta_fd, fasta_file = tempfile.mkstemp()
	fasta_out = open(fasta_file, 'w')
	print >> fasta_out, '>seq\n%s' % seq
	fasta_out.close()

	# print figure to a temp eps file
	eps_fd, eps_file = tempfile.mkstemp()
	weblogo_cmd = 'weblogo --errorbars NO --show-xaxis NO --show-yaxis NO --fineprint "" -c classic -n %d %s < %s > %s' % (len(seq), weblogo_args, fasta_file, eps_file)
        if subprocess.call(weblogo_cmd, shell=True):
            message('Error running weblogo', 'error')

	# copy eps file over and write in my own heights
	start_stack_re = re.compile('^\(\d*\) StartStack')
	out_eps_open = open(out_eps, 'w')
	weblogo_eps_in = open(eps_file)
	line = weblogo_eps_in.readline()
	si = 0
	while line:
		start_stack_match = start_stack_re.search(line)

		# nt column begins
		if start_stack_match:
			print >> out_eps_open, line,

			# loop over 4 nt's
			for i in range(4):
				line = weblogo_eps_in.readline()
				a = line.split()

				nt = a[2][1:-1]
				if nt != seq[si]:
					print >> out_eps_open, line,
				else:
					# change the nt of seq
					a[1] = '%.6f' % heights[si]
					print >> out_eps_open, ' %s' % ' '.join(a)

			# move to next nucleotide
			si += 1
		else:
			print >> out_eps_open, line,

		# advance to next line
		line = weblogo_eps_in.readline()

	# clean
	os.close(fasta_fd)
	os.remove(fasta_file)
	os.close(eps_fd)
	os.remove(eps_file)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
