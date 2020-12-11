'''Common convenience codes for PySCF.'''
from pyscf.scf import ROHF
from pyscf.pbc.scf import KROHF

def load_molcalc(chkfile,mf=ROHF,mfargs=None,dfargs=None):
  ''' Load a cell and MF object from a chkfile.
  Args:
    chkfile: path to chkfile.
    mf: Mean-field object to use. Default is KROHF.
    mfargs: Any internals to mf that aren't in the chkfile, for example, exxdiv, xc, etc.
      They are overwritten by args in chkfile that overlap.
  Returns:
    mol,mf: Cell and mean field object populated by data from chkfile.
  '''
  from pyscf.lib.chkfile import load_mol, load
  if mfargs is None: mfargs = {}
  if dfargs is None: dfargs = {}

  mol = load_mol(chkfile)
  mf = mf(mol,**mfargs).density_fit()
  mf.__dict__.update(load(chkfile,'scf'))
  mf.with_df.__dict__.update(**dfargs)

  return mol,mf

def load_cellcalc(chkfile,mf=KROHF,mfargs=None,dfargs=None):
  ''' Load a cell and MF object from a chkfile.
  Args:
    chkfile: path to chkfile.
    mf: Mean-field object to use. Default is KROHF.
    mfargs: Any internals to mf that aren't in the chkfile, for example, exxdiv, xc, etc.
      They are overwritten by args in chkfile that overlap.
  Returns:
    cell,mf: Cell and mean field object populated by data from chkfile.
  '''
  from pyscf.pbc.lib.chkfile import load_cell, load
  if mfargs is None: mfargs = {}
  if dfargs is None: dfargs = {}

  cell = load_cell(chkfile)
  mf = mf(cell).density_fit()
  mf.__dict__.update(mfargs)
  mf.__dict__.update(load(chkfile,'scf'))
  mf.with_df.__dict__.update(dfargs)

  return cell,mf

