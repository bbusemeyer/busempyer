from pyscf import lib
import sys

assert len(sys.argv)==2,"""
  Usage: chkfile_gaps.py chkfile.
  """

dat=lib.chkfile.load(sys.argv[1],'scf')
nspin=1
print('nspin:',nspin)
if nspin==1:
  egy=dat['mo_energy']
  occ=dat['mo_occ']
  print('Occupations Energies')
  last=0.0
  count=1
  for info in zip(occ,egy):
    o,e=info
    change=e-last
    print("{} {:1.1} {: >07.3f} {: >07.3f}".format(count,o,e,change))
    count+=1
    last=e
else:
  upegy,dnegy=dat['mo_energy']
  upocc,dnocc=dat['mo_occ']
  print('Occupations Energies')
  last=0.0,0.0
  count=1
  for info in zip(upocc,dnocc,upegy,dnegy):
    uo,do,ue,de=info
    change=(ue-last[0],de-last[1])
    print("{} {:1.1} {: >07.3f} {: >07.3f}".format(count,uo,ue,change[0]))
    count+=1
    last=ue,de
