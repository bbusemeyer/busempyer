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
    queue="batch",
    name="%s/%s.out"%(cwd,base),
    time="72:00:00",
    nn=1,np=8,
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
  #"lambda":[0.01,0.02,0.03,0.1,0.2]
  #"kpoint":[2,4,6,8,10],
  #"ecutwfc":[50,60,70,80,90,100],
  #"ecutrho":[400]
  #"psuedo":[
  #  ("Si","Si.pbe-n-rrkjus_psl.0.1.UPF"),
  #  ("Si","Si.rel-pbe-n-rrkjus_psl.0.1.UPF"),
  #  ("Si","Si.pz-vbc.UPF")
  "nqx":[1,2,3]
}

# Now's the part where you mess with something.
for key in changes.keys():
  lines = deepcopy(baselines)
  if key=="kpoint": 
    for newval in changes["kpoint"]:
      lines[-1] = [str(newval),str(newval),str(newval),"0","0","0"]
      base = "conv_%s_%s"%(key,newval)
      submit_job(base,lines,cwd)
    continue
  if key=="psuedo":
    start,end = 0,0
    # Find where pseudos are chosen.
    for li, line in enumerate(lines):
      if "ATOMIC_SPECIES" in line:
        start = li+1
    for li, line in enumerate(lines):
      if "ATOMIC_POSITIONS" in line:
        end = li-1
    # Replace for every species.
    for atom,pot in changes[key]:
      poss = []
      for li in range(start,end):
        if atom in lines[li][0]:
          poss.append(li)
      for pos in poss:
        lines[pos][-1] = pot
        base = "pseudo_%s"%pot
        submit_job(base,lines,cwd)
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
  if key=="nqx":
    for newval in changes[key]:
      for line in lines:
        if any(["nqx" in word for word in line]):
          line[-1] = str(newval)
      base = "conv_%s_%s"%(key,newval)
      submit_job(base,lines,cwd)

  # Basic key replacement.
  for line in lines:
    if key in line:
      for newval in changes[key]:
        line[-1] = str(newval)
        base = "conv_%s_%s"%(key,newval)
        submit_job(base,lines,cwd)
