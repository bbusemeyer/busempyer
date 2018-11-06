#!/usr/bin/python
import numpy as np

def safeint(thing):
  try: return int(thing)
  except ValueError: return thing

def read_map(inpf):
  ''' Generate a nearest neighbor table from a crystal output file.

  Args:
    inpf (file object): Crystal output file opened for reading.
  Returns:
    dict: Map of ints (atom numbers labeled by crystal, which is 1-based indexing) to a list
      [neighbor_number, neighbor_element, ix,iy,iz, distance].
      The ix etc. label if it crosses a periodic bounary.
      The information is similar to how it appears in the crystal output.
  '''
  line = "Start"
  found_spot = False
  nnmap = {}
  while line!='':
    line = inpf.readline()
    if not found_spot:
      if "NEIGHBORS (ATOM LABELS AND CELL INDICES)" in line:
        found_spot = True
      continue
    if line == '\n': continue
    if "SYMMETRY" in line: break
      
    anum, elem, nnbr, rnbr, dumy = line.split()[:5]
    anum=int(anum)
    nnbr=int(nnbr)

    tmp = line[33:-1]
    for nlines in range(nnbr//3 - 1):
      tmp += inpf.readline()[:-1]
    if (nnbr > 3) and (nnbr % 3 != 0):
      tmp += inpf.readline()[:-1]
    tmp = tmp.replace('-',' -').split()
    nbrs = [[safeint(t) for t in tmp[5*i:5*(i+1)]] + [float(rnbr)] for i in range(nnbr)]

    if anum in nnmap.keys():
      nnmap[anum] += nbrs
    else:
      nnmap[anum] = nbrs
  return nnmap

def map_to_table(nnmap):
  ''' Convert the mapping (good at looking up neighbors of an atom) to a table (good at looking at all distances).'''
  nntable = []
  for key in nnmap.keys():
    for neighbors in nnmap[key]:
      nntable.append([key]+neighbors)
  return nntable

if __name__ == "__main__":
  import sys
  inpf = open(sys.argv[1],'r')
  #outf = open("neighbor_table.csv",'w')
  
  nnmap = read_map(inpf)
  print(nnmap)

  ## Output data to usable table.
  #print("anum  nnum  nelm  cx  cy  cz  nrad\n")
  #for key in nnmap.keys():
  #  for neighbors in nnmap[key]:
  #    print('%s  %s  %s  %s  %s  %s  %s\n'%tuple([key]+neighbors))
