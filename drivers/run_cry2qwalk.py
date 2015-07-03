#!/usr/bin/python

from cryfiles_io import gen_cry2qwalk
from subprocess import check_output
import sys

for dftfn in sys.argv[1:]:
  gen_cry2qwalk(dftfn)
  dr = '/'.join(['.']+dftfn.split('/')[:-1])
  print check_output('qsub {0}/qsub.in'.format(dr),shell=True)
