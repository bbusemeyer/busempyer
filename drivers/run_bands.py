import cryfiles_io as cio
import sys
import subprocess as sub

# I'm not sure how general these are, but they're for FeSe.
cryinp = open(sys.argv[1],'r')
natoms = [4,4]
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
with open("band.inp",'w') as bandf:
  bandf.write(cio.gen_properties(cryinp,natoms,kpath,denom,projs,above=10,below=10))
#print sub.check_output("properties < band.inp > band.out",shell=True)
