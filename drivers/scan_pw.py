import sys
import os
from copy import deepcopy
import subprocess as sp
from mython import gen_qsub

def writelines(lines):
  strlines = [" ".join(line) for line in lines]

def submit_job(base,lines,cwd):
  if not os.path.exists(base): os.mkdir(base)
  os.chdir(base)
  with open("%s.inp"%base,'w') as outf:
    outf.write("\n".join([" ".join(l) for l in lines]))
  pc = ["module load openmpi/1.4-gcc+ifort"]
  qin = gen_qsub("~/bin/pw.x < %s.inp"%(base),
    stdout="%s.out"%(base),
    queue="physics",
    name="%s/%s.out"%(cwd,base),
    time="72:00:00",
    nn=1,
    prep_commands=pc)
  print sp.check_output("qsub %s"%qin,shell=True)
  os.chdir(cwd)

if len(sys.argv) < 2:
  print "You need to enter a base file, dummy!"
  exit(1)
else:
  basefn = sys.argv[1]

baselines = []
with open(basefn,'r') as basef:
  for line in basef:
    baselines.append(line.split())
# Make base file for easy comparison.
with open("base.inp",'w') as outf:
  outf.write("\n".join([" ".join(l) for l in baselines]))

cwd = os.getcwd()

changes = {
    "lambda":[0.01,0.02,0.03,0.1,0.2]
  #"kpoint":[2,4,6,8,10],
  #"ecutwfc":[50,60,70,80,90,100],
  #"ecutrho":[400]
}

# Now's the part where you mess with something.
for key in changes.keys():
  lines = deepcopy(baselines)
  if key=="kpoint": 
    for newval in changes[key]:
      lines[-1] = [str(newval),str(newval),str(newval),"0","0","0"]
      base = "conv_%s_%s"%(key,newval)
      submit_job(base,lines,cwd)
    continue
  if key=="ecutwfc":
    ecutwfc,ecutrho=0,0
    for li,line in enumerate(lines):
      if "ecutwfc" in line: ecutwfc=li
      if "ecutrho" in line: ecutrho=li
    for newval in changes[key]:
      lines[ecutwfc][-1] = str(newval)
      lines[ecutrho][-1] = str(10*newval)
      base = "conv_%s_%s"%(key,newval)
      submit_job(base,lines,cwd)
    continue # TODO: test new configuration.
  # Basic key replacement.
  for line in lines:
    if key in line:
      for newval in changes[key]:
        line[-1] = str(newval)
        base = "conv_%s_%s"%(key,newval)
        submit_job(base,lines,cwd)
