'''Common convenience codes for PySCF.'''
from pyscf.pbc.scf import KRHF

def load_cellcalc(chkfile,mf=KRHF,mfargs=None,dfargs=None):
  ''' Load a cell and MF object from a chkfile.
  Args:
    chkfile: path to chkfile.
    mf: Mean-field object to use. Default is KRHF.
    mfargs: Any internals to mf that aren't in the chkfile, for example, exxdiv, xc, etc.
      They are overwritten by args in chkfile that overlap.
  Returns:
    cell,mf: Cell and mean field object populated by data from chkfile.
  '''
  from pyscf.pbc.lib.chkfile import load_cell, load
  if mfargs is None: mfargs = {}
  if dfargs is None: mfargs = {}

  cell = load_cell(chkfile)
  mf = mf(cell,**mfargs).density_fit()
  mf.__dict__.update(load(chkfile,'scf'))
  mf.with_df.__dict__.update(**dfargs)

  return cell,mf

