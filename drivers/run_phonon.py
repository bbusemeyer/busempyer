import qefiles_io as qeio
import os
import mython as my
import sys
import subprocess as sub

qpts = [(0.0,0.0,0.0),(0.25,0.0,0.0),(0.5,0.0,0.0),(0.75,0.0,0.0)]
scfinp = open(sys.argv[1],'r')
inps = qeio.gen_phonon(qpts,scfinp,bundle=True)

cwd = os.getcwd()
for fname in inps:
  if '/' in fname:
    dname = '/'.join(fname.split('/')[:-1])
    fname = fname.split('/')[-1]
    os.chdir(dname)

  qsub = my.gen_qsub(
    '~/bin/ph.x < %s'%fname,
    stdout=fname.replace('inp','out'),
    loc=os.getcwd(),
    nn=1,np=8)
  print sub.check_output("qsub %s"%qsub,shell=True)
  os.chdir(cwd)
