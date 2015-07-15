from qfiles_io import gen_dmc
from subprocess import check_output
from numpy import array
import sys

# Real k-points
kpoints = array([1,4,17,20,93,96,109,112]) - 1

dftfn = sys.argv[1]
sysfns = [dftfn.replace('.d12','_%d.sys'%k) for k in kpoints]

time=60*2
nproc=512*4

qsub = []
qsub.append('qsub')
qsub.append('-q prod')
qsub.append('-A SuperMatSim')
qsub.append('-t {time}'.format(time=time))
qsub.append('-n {nproc}'.format(nproc=nproc))
qsub.append('--mode c32')
qsub.append('-o {gamma}.dmc.out'.format(gamma='/'.join((loc,gamma))))
qsub.append('~/bin/qwalk')
for loc,root in info:
  qsub.append(loc+'/'+root+'.dmc')
qin = ' '.join(qsub)

qin = gen_dmc(sysfns)
check_output(qin,shell=True)
#print qin
