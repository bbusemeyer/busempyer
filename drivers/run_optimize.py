#!/usr/bin/python

from cryfiles_io import gen_cry2qwalk
from qfiles_io   import gen_optimize
from subprocess import check_output
import sys

NRESARTS = 4

for dftfn in sys.argv[1:]:

  qin = gen_optimize(dftfn)
  job = check_output('qsub {0}'.format(qin),shell=True)
  for redo in range(NRESARTS-1):
    job = check_output('qsub {0} -W depend=afterok:{1}'.format(qin,job),
                       shell=True)
    print job
