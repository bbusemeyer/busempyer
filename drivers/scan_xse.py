#!/usr/bin/python
from mython import lines2str,gen_qsub
from subprocess import check_output
import os
import sys

if len(sys.argv) != 2:
  print "Wrong number of arguments, buddy."
  exit()
else:
  print "Basing input off of {base}".format(base=sys.argv[1])

#dx_list = [-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0.03,0.02,0.01]
dx_list = [-0.06,-0.05,-0.04,-0.03,0.03,0.02]
hyb_list = [0]#,2,4,6,8,10,15,20]

baselines = []
with open(sys.argv[1],'r') as basef:
  for line in basef:
    baselines.append(line.split())

spot = (None,None)
val  = None
for li in range(len(baselines)):
  line = baselines[li]
  if "END" == line[0]:
    print "Didn't find spot."
    break
  else:
    if "234" == line[0]:
      spot = (li,3)
      val  = float(line[3])
      break
for li in range(len(baselines)):
  line = baselines[li]
  if [] == line:
    print "Didn't find hybrid spot."
    break
  else:
    if "HYBRID" == line[0]:
      hyb_spot = li+1
      break

with open('base.d12','w') as outf:
  outf.write(lines2str(baselines))

for hy in hyb_list:
  for dx in dx_list:
    dr = 'hy_%d_dx_%.2f'%(hy,dx)
    baselines[spot[0]][spot[1]] = val+dx
    baselines[hyb_spot][0] = hy
    if not os.path.exists(dr):
      os.makedirs(dr)
    with open('{dr}/{dr}.d12'.format(dr=dr),'w') as outf:
      outf.write(lines2str(baselines))

    loc = '/'.join((os.getcwd(),dr))
    pc  = ['module load openmpi/1.4-gcc+ifort',
           'cp {dr}.d12 INPUT'.format(dr=dr)]
    fc = ['rm *.pe[0-9] *.pe[0-9][0-9]']
    qin = gen_qsub('~/bin/Pcrystal',
                   loc=loc,
                   name='_'.join((sys.argv[1][:3],dr)),
                   stdout=dr+'.d12.out',
                   time='4:00:00',
                   nn=2,np=12,
                   queue='secondary',
                   prep_commands = pc,
                   final_commands = fc)
    #print check_output('cp -v {0}/fort.79 {0}/fort.20'.format(dr),shell=True)
    print check_output('cp -v fort.9 {0}/fort.20'.format(dr),shell=True)
    print check_output('qsub {0}'.format(qin),shell=True)
