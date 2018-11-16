import numpy as np

def read_xyz(xyzfn):
  lines=open(xyzfn,'r').read().split('\n')
  index={}
  atoms={}
  for line in lines[2:]:
    if line=='':continue
    words=line.split()
    atom=words[0]
    if atom in index.keys():
      index[atom]+=1
    else:
      index[atom]=1
    atoms[(atom,index[atom])]=np.array(words[1:],dtype=float)
  return atoms

def export_xyz(atoms,outfn='new.xyz'):
  outlines = [
    str(len(atoms.keys())),
    "Made by patch_structure.py",
  ]
  for key in sorted(atoms.keys()):
    outlines += [
        "{:4} {: 10.10f} {: 10.10f} {: 10.10f}".format(key[0],*atoms[key])
    ]
  with open(outfn,'w') as outf:
    outf.write('\n'.join(outlines))
