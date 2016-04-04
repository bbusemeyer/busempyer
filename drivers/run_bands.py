import cryfiles_io as cio
import sys
import subprocess as sub

# I'm not sure how general these are, but they're for FeSe.
cryinp = open(sys.argv[1],'r')
natoms = [2,2]
kpath = [
    [0,0,0],
    [1,0,0],
    [1,1,0],
    [1,1,1],
    [1,0,1],
    [0,0,1]
  ]
denom = 2
projs = []
bandout = open("band.inp",'w')
bandout.write(cio.gen_properties(cryinp,natoms,kpath,denom,projs))
bandout.close()
print sub.check_output("properties < band.inp > band.out",shell=True)
