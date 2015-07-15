from qfiles_io import gen_dmc
from subprocess import check_output
from numpy import array
import sys

# Real k-points
kpoints = array([1,4,17,20,93,96,109,112]) - 1

dftfn = sys.argv[1]
sysfns = [dftfn.replace('.d12','_%d.sys'%k) for k in kpoints]

qins = gen_dmc(sysfns,'~/bin/gosling')
for qin in qins:
  print check_output(qin,shell=True)
