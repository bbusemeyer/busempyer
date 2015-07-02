#!/usr/bin/python
import sys
from numpy import *

inpf = open(sys.argv[1],'r')
outf = open("neighbor_table.csv",'w')
line = "Start"
found_spot = False

# Read and organize data.
nntable = {}
while line!='':
  line = inpf.readline()
  if not found_spot:
    if "NEIGHBORS (ATOM LABELS AND CELL INDICES)" in line:
      found_spot = True
    continue
  if line == '\n': continue
  if "SYMMETRY" in line: break
    
  anum, elem, nnbr, rnbr, dumy = line.split()[:5]
  anum, nnbr = map(int,[anum,nnbr])
  #rnbr = float(rnbr)
  tmp = line[33:-1]
  for nlines in range(nnbr//3 - 1):
    tmp += inpf.readline()[:-1]
  if (nnbr > 3) and (nnbr % 3 != 0):
    tmp += inpf.readline()[:-1]
  tmp = tmp.replace('-',' -').split()
  nbrs = [tmp[5*i:5*(i+1)] + [rnbr] for i in range(nnbr)]

  if anum in nntable.keys():
    nntable[anum] += nbrs
  else:
    nntable[anum] = nbrs

# Output data to usable table.
outf.write("anum,nnum,nelm,cx,cy,cz,nrad\n")
for key in nntable.keys():
  for neighbors in nntable[key]:
    outf.write('%d,%s,%s,%s,%s,%s,%s\n'%tuple([key]+neighbors))
inpf.close()
outf.close()
