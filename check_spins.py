import pandas as pd
import numpy as np
from pyscf import lib
from pyscf.scf.uhf import mulliken_meta
from pyscf.dft.uks import UKS

mol=lib.chkfile.load_mol('pyscf_driver.py.chkfile')
scf=UKS(mol)
dm=scf.from_chk('pyscf_driver.py.chkfile')

orbdf=pd.DataFrame(mol.spherical_labels(),columns=['atid','at','orb','type'])

pops=mulliken_meta(mol,dm,verbose=1)
orbdf['up']=pops[0][0]
orbdf['dn']=pops[0][1]
orbdf['spin']=orbdf['up']-orbdf['dn']
orbdf['charge']=orbdf['up']+orbdf['dn']

atdf=orbdf.groupby(['atid','at']).agg({
    'spin':np.sum,
    'charge':np.sum,
    'up':np.sum,
    'dn':np.sum
  })

print(atdf[['up','dn','charge','spin']])

print("Max spin:",abs(atdf['spin']).max())
print("Total charge:",atdf['charge'].sum())
print("Total spin:",atdf['spin'].sum())
