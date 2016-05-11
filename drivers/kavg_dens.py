# Brian Busemeyer 04/2015
# Extract up and down orbitals
# Compute twist-averaged density.
# Output as cube file.
# Twist average over set of dmc results input.

from numpy import array
import multiprocessing as mp
import cubetools as ct
import sys
from operator import add
from imp import reload
reload(ct)

if len(sys.argv) < 3:
  print("Arguements : [list of files to average] outputfile")
  raise AssertionError("Not enough arguments.")

flist = sys.argv[1:-1]
w_k = [1./len(flist) for i in flist] # Assumes equal k-point weight
print("weights:", w_k)

def reader(fi): return ct.read_cube(open(fi,'r'),qwalk_patch=True)
with mp.Pool(8) as pool:
  dens_k = pool.map(reader,flist)

print("Output to "+sys.argv[-1])
dens = dens_k[0]
dens['comment'] = "Twist-average density\n"
dens['data'] = w_k[0]*dens_k[0]['data']
for i in range(1,len(w_k)):
  dens['data'] += w_k[i]*dens_k[i]['data']
ct.write_cube(dens,open(sys.argv[-1],'w'))
