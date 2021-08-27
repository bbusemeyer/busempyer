''' Some drivers for ASE (atomic simulation environment) tools.

Some useful ASE routines:
  - Atoms.get_all_distances(mic=True): Get distance table, respecting periodic boundary conditions.

'''

import ase
import numpy as np


def ase_from_crystal(crystaloutfn):
  ''' Create an ASE Atoms object by reading structure from a crystal output file.

  Args:
    crystaloutfn (str): output file name.
  Returns:
    Atoms: structure this represents.
  '''
  # Assumes that crystal puts all atoms of the same element togther.
  # If this causes bugs probably just need an argsort somewhere.
  reader=open(crystaloutfn,'r')
  lines=reader.readlines()

  lparms={}
  pos={'elemnums':[],'coords':[]}

  for lnum,line in enumerate(lines):
    if 'LATTICE PARAMETERS' in line and 'BOHR' in line:
      # Sometimes its in two lines, sometimes its in three.
      try:
        lparms=dict(zip(('a','b','c','alpha','beta','gamma','vol'),(float(n) for n in lines[lnum+2].split())))
      except ValueError:
        lparms=dict(zip(('a','b','c','alpha','beta','gamma','vol'),(float(n) for n in lines[lnum+3].split())))
    if 'CARTESIAN COORDINATES' in line:
      cursor=lnum+4
      while True:
        try:
          atnum,elemnum,elem,x,y,z = lines[cursor].split()
          elemnum=int(elemnum)-200 # Assumes pseudopot.
          x,y,z = float(x),float(y),float(z)
          pos['elemnums'].append(elemnum)
          pos['coords'].append((x,y,z))
          cursor+=1
        except ValueError:
          break
    # Once we're done, stop reading because the file can be huge.
    if (len(lparms)*len(pos['elemnums'])*len(['coords']))>0:
      break

  return ase.Atoms(pos['elemnums'],pos['coords'],pbc=True,cell=[lparms[k] for k in ('a','b','c','alpha','beta','gamma')])

def get_nnmap(atoms):
  ''' Make a map of atom indicies to neighbors, ordered by closeness.
  
  Args:
    atoms (Atoms): ASE Atoms object.
  Returns:
    dict: nnmap[atomnum] = ( (atom number, distance) for each other atom ) in order of decreasing closeness.
  '''
  distances=atoms.get_all_distances(mic=True)
  nnmap={}
  for atomnum,atomdist in enumerate(distances):
    order=atomdist.argsort()
    nnmap[atomnum]=list(zip(
        np.arange(atomdist.shape[0])[order][1:],
        np.array(atoms.get_chemical_symbols())[order][1:],
        atomdist[order][1:]
      ))
  return nnmap
