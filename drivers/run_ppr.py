from qfiles_io import gen_ppr
from subprocess import check_output
from numpy import array
from mython import gen_qsub
import sys
import os

if len(sys.argv) < 2:
  print "Input is DFT file names that generates the ppr input." 
  exit()

# Real k-points
kpoints = array([1]) - 1 # 2x2x2
#kpoints = array([1,3,8,10,27,29,34,36]) - 1 # 2x2x2
#kpoints = array([1,4,17,20,93,96,109,112]) - 1 # 4x4x4
#kpoints = array([1,5,30,34,227,231,256,260]) - 1

wd = os.getcwd()
for dftloc in sys.argv[1:]:
  dftfn  = dftloc.split('/')[-1]
  roots  = [dftfn.replace('.d12','_{k}'.format(k=k)) for k in kpoints]
  gamma  = roots[0]
  loc    = '/'.join(dftloc.split('/')[:-1])

  os.chdir(wd + '/' + loc)

  pprfns = [gen_ppr(root,gosling="/home/busemey2/bin/gosling",
                         jast=gamma+'.opt.jast2') for root in roots]
    
  pc = []
  pc += ['module load openmpi/1.4-gcc+ifort']
  pc += ['module load gcc/4.7.1']
  pc += ['module load intel/11.1']
  for pprfn in pprfns:
    qin = gen_qsub('~/bin/qwalk {0}'.format(pprfn),
                   stdout=pprfn+'.out',
                   name='/'.join((os.getcwd(),pprfn)),
                   time='00:01:00',
                   nn=2,np=12,
                   queue='test',
                   prep_commands=pc)
    print check_output('qsub '+qin,shell=True)
    #print qin

  os.chdir(wd)
