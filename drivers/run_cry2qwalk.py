#!/usr/bin/python

from cryfiles_io import gen_cry2qwalk
from subprocess import check_output
import sys

for dftfn in sys.argv[1:]:
  qin = gen_cry2qwalk(dftfn)
  print check_output('qsub {0}'.format(qin),shell=True)
