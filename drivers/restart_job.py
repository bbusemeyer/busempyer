
import mython as my
import cryfiles_io as cio
import shutil as sh
import subprocess as sub
import os
import sys

# qsub options.
exe = "~/bin/Pcrystal" 
nn  = 2
time = "04:00:00"
queue = "secondary"
pc = ["module load openmpi/1.4-gcc+ifort","rm INPUT","cp XXX INPUT"]
fc = ["rm *.pe[0-9]","rm *.pe[0-9][0-9]"]

use_fort = "../fort.9"

inpfns = open(sys.argv[1],'r').readlines()
inpfs = [open(i[:-1],'r') for i in inpfns]

cwd = os.getcwd()

for inpf in inpfs:
  loc = "/".join(inpf.name.split("/")[:-1])
  inpfn = inpf.name.split("/")[-1]
  needs_guessp = True
  os.chdir(loc)

  inpstr = ""
  for line in inpf:
    inpstr += line
    if "GUESSP" in line: needs_guessp=False
  inplines = inpstr.split("\n")

  if needs_guessp:
    i=0
    while inplines[-i] != "END":
      i += 1
    inplines[-i] = "GUESSP"
    try:
      inplines[-i+1] = "END"
    except IndexError:
      inplines.append("END")

  newinpfn = "new."+inpfn
  with open(newinpfn,'w') as outf:
    outf.write("\n".join(inplines))

  sh.copyfile(use_fort,"fort.20")
  pc[-1] = "cp %s INPUT"%newinpfn
  qsub = my.gen_qsub(exe,
    stdout = newinpfn+".out",
    loc = os.getcwd(),
    name = "%s %s"%(os.getcwd(),newinpfn),
    time = time,
    nn = nn,
    queue = queue,
    prep_commands = pc,
    final_commands = fc
  )
  print newinpfn, sub.check_output("qsub %s"%qsub,shell=True)

  os.chdir(cwd)
