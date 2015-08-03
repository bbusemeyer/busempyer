from qfiles_io import gen_dmc
from subprocess import check_output
from numpy import array
from mython import gen_qsub
import sys
import os

system = ['mira','taub'][1]

if len(sys.argv) != 2:
  print "Input is DFT file name that generates the DMC input." 
  exit()

# Real k-points
kpoints = array([1,4,17,20,93,96,109,112]) - 1

dftloc = sys.argv[1]
dftfn  = dftloc.split('/')[-1]
roots  = [dftfn.replace('.d12','_{k}'.format(k=k)) for k in kpoints]
gamma  = roots[0]
loc    = '/'.join(dftloc.split('/')[:-1])

wd = os.getcwd()
os.chdir(wd + '/' + loc)

dmcfns = [gen_dmc(root,jast=gamma+'.opt.jast2') for root in roots]
  
# qsub routine for Mira.
if system == 'mira':
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
  qsub += dmcfns
  qin = ' '.join(qsub)
  #print check_output(qin,shell=True)
  print qin
elif system == 'taub':
  pc =  ['module load openmpi/1.4-gcc+ifort']
  for dmcfn in dmcfns:
    qin = gen_qsub('~/bin/qwalk {0}'.format(dmcfn),
                   stdout=dmcfn+'.out',
                   name='/'.join((os.getcwd(),dmcfn)),
                   time='4:00:00',
                   nn=4,np=12,
                   queue='secondary',
                   prep_commands=pc)
    print check_output('qsub '+qin,shell=True)
    #print qin
else:
  print "You need to add appropriate submission routine."

os.chdir(wd)
