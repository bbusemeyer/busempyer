''' Some standard routines for PySCF IO and queue interaction.'''
import os, yaml, shutil, ccq_sub_py
from pyscf.gto import Mole
from pyscf.pbc.gto import Cell
from pyscf.pbc.lib.chkfile import save_cell
from pyscf.lib.chkfile import load_mol, save_mol
from h5py import File
from pyscf.pbc.tools.pyscf_ase import ase_atoms_to_pyscf

# Filename convention. f"{SCFNAME}.py" will be run.
SCFNAME = "scfcalc"

def make_cell(atomase,**cellargs):
  ''' Make a cell from an ASE Atoms object.'''

  cell = Cell()
  cell.atom = ase_atoms_to_pyscf(atomase)
  cell.a = atomase.cell.tolist()
  cell.build(**cellargs)

  return cell

def runcalc(loc,cell,mfargs={},qsubargs={'time':'6:00:00','queue':'gen'},guess=None,dfints=None,meta={}): 
  ''' Deposit run input into a location and run.
  Args:
    loc (str): directory for pyscf.
    cell (Cell): PySCF Mole or Cell.
    mfargs (dict): non-default args for mean field object.
  Returns:
    bool: if the directory was newly prepped.
  '''
  print(f"\n --- Checking job {loc}. ---")
  if loc[-1] != '/': loc+='/'
  cwd = os.getcwd()

  if os.path.exists(f"{loc}{SCFNAME}.py") or os.path.exists(f"{loc}{SCFNAME}.yaml"):
    print("Already started.")
    #yaml.dump(meta,stream=open(f"{loc}meta.yaml",'w'))
    return False

  if not os.path.exists(loc): os.mkdir(loc)

  print(f"Preparing calculation in {loc}{SCFNAME}...")

  if type(cell) == Cell:
    save_cell(cell,f"{loc}{SCFNAME}.chk")
  elif type(cell) == Mole:
    save_mol(cell,f"{loc}{SCFNAME}.chk")
  else:
    raise AssertionError("Struture type not recognized.")

  yaml.dump(mfargs,stream=open(f"{loc}{SCFNAME}.yaml",'w'))
  yaml.dump(meta,stream=open(f"{loc}meta.yaml",'w'))

  if dfints is not None:
    shutil.copyfile(dfints,f"{loc}{SCFNAME}_gdf.h5")

  if guess is not None:
    shutil.copyfile(guess,f"{loc}guess.chk")

  shutil.copyfile(f"{SCFNAME}.py",f"{loc}{SCFNAME}.py")

  os.chdir(loc)
  ccq_sub_py.qsub(SCFNAME+'.py',**qsubargs)
  os.chdir(cwd)

  print(f"Done running {loc}.")
  return True

def readcalc(loc):
  if loc[-1] != '/': loc+='/'
  print(f"\nReading results from {loc}{SCFNAME}.")
  results = {}
  results['loc'] = loc

  root = loc+SCFNAME
  scfyaml = root+'.yaml'
  scfchk  = root+'.chk'
  meta    = loc+'meta.yaml'
  stdout  = root+'.py.out'

  if os.path.exists(scfyaml):
    results.update(yaml.safe_load(open(scfyaml,'r')))

  if os.path.exists(meta):
    results['meta'] = yaml.safe_load(open(meta,'r'))

  if os.path.exists(scfyaml):
    results.update(yaml.safe_load(open(scfyaml,'r')))

  if os.path.exists(scfchk):
    struct = load_mol(scfchk)
    results['basis'] = struct.basis
    results['ecp'] = struct.ecp
    results['spin'] = struct.spin
    if 'exp_to_discard' in struct.__dict__: results['exp_to_discard'] = struct.exp_to_discard
    else:                                   results['exp_to_discard'] = None
    results['atoms'] = [struct.atom_symbol(i) for i in range(struct.natm)]
    results['nelec'] = struct.tot_electrons()
    results['verbose'] = struct.verbose

    #output = pyscf.pbc.lib.chkfile.load(scfchk,'scf')
    output = File(scfchk)
    if 'scf' in output.keys():
      results['e_tot'] = output['scf']['e_tot'][()]
      print("Found energy:",results['e_tot'])
    else: print("No results found.")

  results['converged'] = False
  if os.path.exists(stdout) and 'converged SCF energy' in open(stdout).read():
    results['converged'] = True

  return results
