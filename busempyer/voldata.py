''' Collection of tools for reading and formatting 3d data.'''
from numpy import array, asarray, diag, zeros, empty, dot, eye, cross, meshgrid
from numpy.linalg import inv
from pyscf.pbc.dft.numint import eval_rho

def main():
  ''' Defines a command-line interface.'''
  import argparse
  parser=argparse.ArgumentParser("data3d")
  #parser.add_argument('inpfn',type=str,
  #    help='Input to be submitted.')
  #parser.add_argument('-t',dest='time',default='1:00:00',type=str,
  #    help='Time string.'+dft)
  #parser.add_argument('-q',dest='queue',default='ccq',type=str,
  #    help='Queue.'+dft)
  #parser.add_argument('-l',dest='local',action='store_true',
  #    help='Run locally instead of on the cluster'+dft)
  args=parser.parse_args()

  print(args)

# TODO decide units.
class VolData:
  def __init__(self,data=None,grid=None,voxel=None,positions=None,latvecs=None,origin=(0.0,0.0,0.0),meta=None,verbose=False):
    ''' Contains data and routines for exporting data for volumetric data formats.
    Args:
      data (array-like): 3-D data in a list.
      grid (array-like): Coordinates in 3D space of the points. 
      voxel (array-like): 3x3 array defines edges of the voxel.
      positions (list): [('element',x,y,z) for atom in atoms].
      latvecs (array-like): 3x3 array defines the periodicity.
      origin (array-like): Origin of data appears here in plotted data.
      meta (dict): Any special information to store.
      verbose (bool): Print progress for larger jobs.
    '''
    self.data       = data if data is not None else zeros((3,3,3))
    self.npoints    = self.data.shape
    self.grid       = grid
    self.voxel      = voxel if voxel is not None else diag(array(self.data.shape,dtype=float)**-1)
    self.positions  = positions if positions is not None else []
    self.latvecs    = latvecs
    self.origin     = asarray(origin)
    self.meta       = meta if meta is not None else {}
    self.verbose    = verbose

    for atom in self.positions: 
      atom[1] = asarray(atom[1])

  # Can add from_cube from cubetools easily.

  def load_from_cell_(self, cell, npoints, grid=None, origin=(0.0,0.0,0.0)):
    ''' Load cell geometry data, and make grid in that cell.'''
    self.npoints = asarray(npoints)
    self.latvecs = asarray(cell.lattice_vectors())
    self.voxel = self.latvecs/self.npoints[:,None]
    self.origin = asarray(origin)
    self.grid, fracoords = make_grid(self.voxel,self.npoints,origin=origin)

    cart = array([cell.atom_coord(i) for i in range(cell.natm)]).T
    if sum(self.origin) > 1e-10:
      cart = _set_nonzero_origin(cart, self.latvecs, self.origin)
    self.positions = [(cell.atom_symbol(i),cart[:,i]) for i in range(cell.natm)]

  def compute_pyscf_points_(self, mol, coeff_or_dm, dtype='orbital', verbose=False):
    self.data = compute_pyscf_points(mol, self.grid, coeff_or_dm, dtype, verbose).reshape(self.npoints)

  def from_pyscf(self, cell, coeff_or_dm, npoints, dtype='orbital', origin=(0.0,0.0,0.0)):
    ''' Convenience function to load and compute all data from a PySCF run.'''

    self.load_from_cell_(cell, npoints, origin)
    self.data = compute_pyscf_points(cell, self.grid, coeff_or_dm, dtype, verbose=self.verbose).reshape(npoints)

    return self

  # See XSF format specs on "general grids" for extra printing of 0-index elements.
  def write_xsf(self,outf):
    from ase.units import Bohr
    ''' Write internal data into an XSF format.
    Note: units of XSF are Angstrom!
    Currently requires periodic systems, but can easily be generalized to non-PBC
    '''
    if self.verbose: print("\nExporting XSF.")
    if type(outf)==str: outf=open(outf,'w')
    outf.write("CRYSTAL\n")
    outf.write("PRIMVEC\n")
    natoms=len(self.positions)
    for i in range(0,3):
      outf.write("  {0: 20.16e} {1: 20.16e} {2: 20.16e}\n".format(*(self.latvecs[i,:]*Bohr)))
    outf.write("PRIMCOORD\n")
    outf.write("  %i 1\n"%natoms)
    for i in range(0,natoms):
      outf.write("  {sym} {0: 20.16e} {1: 20.16e} {2: 20.16e}\n"\
          .format(*(self.positions[i][1]*Bohr),sym=self.positions[i][0]))
    outf.write("BEGIN_BLOCK_DATAGRID_3D\n  from_VolData \n")
    outf.write("BEGIN_DATAGRID_3D\n")
    outf.write("  {0:10} {1:10} {2:10}\n".format(*[s+1 for s in self.data.shape]))
    outf.write("  0.0 0.0 0.0\n")
    for i in range(0,3):
      outf.write("  {0: 20.16e} {1:20.16e} {2:20.16e}\n".format(*(self.data.shape[i]*self.voxel[i,:]*Bohr)))
    
    count=0
    outf.write('  ')
    # Note see "general grid" in XSF definition format.
    # Essentially, need to copy the 0 component again after finishing each loop.
    for z in range(0,self.data.shape[2]):
      for y in range(0,self.data.shape[1]):
        for x in range(0,self.data.shape[0]):
          count = _write3d(self.data[x,y,z],outf,count)
        count = _write3d(self.data[0,y,z],outf,count)
      for x in range(self.data.shape[0]):
        count = _write3d(self.data[x,0,z],outf,count)
      count = _write3d(self.data[0,0,z],outf,count)
    for y in range(0,self.data.shape[1]):
      for x in range(self.data.shape[0]):
        count = _write3d(self.data[x,y,0],outf,count)
      count = _write3d(self.data[0,y,0],outf,count)
    for x in range(self.data.shape[0]):
      count = _write3d(self.data[x,0,0],outf,count)
    count = _write3d(self.data[0,0,0],outf,count)
    outf.write('\n')
    outf.write("END_DATAGRID_3D\n")
    outf.write("END_BLOCK_DATAGRID_3D\n")

  def integrate(self, gridsel=None, transform=None):
    ''' Integrate the data by using a piecewise constant approximation.
    Args:
      gridsel: boolean array of same shape as self.grid for selecting the points for plotting.
      transform: call this on the data before integrating, e.g. abs(). 
      '''
    voxel_volume = cross(self.voxel[0], self.voxel[1]) @ self.voxel[2]

    data = self.data \
        if gridsel is None else self.data.ravel()[gridsel]
     
    return data.sum()*voxel_volume \
        if transform is None else transform(data).sum()*voxel_volume

  def subtract_(self, other):
    ''' Subtract two VolData object's data and store in this data.'''
    self.data -= other.data

def _write3d(data,outf,count):
  ''' data format for writing volume data.'''
  outf.write("{0: 20.16e} ".format(data))
  count+=1
  if count%5==0:
    outf.write('\n  ')

  return count

def compute_pyscf_points(mol, grid, data, dtype='orbital', verbose=False):
  ''' Compare orbital values at points npoints multiples of latvecs'''
  from pyscf import lib
  from pyscf.pbc.gto import Cell

  GTOval = 'GTOval'
  if isinstance(mol, Cell):
    GTOval = 'PBC' + GTOval

  # Compute density on the .cube grid
  ngrids = grid.shape[0]

  blksize = min(8000, ngrids)
  data_on_grid = empty(ngrids)

  if 'orb' in dtype:
    for ip0, ip1 in lib.prange(0, ngrids, blksize):
      if verbose: print(f"Done with {ip0}/{ngrids} = {ip0/ngrids:0.2f}.")
      ao = mol.eval_gto(GTOval, grid[ip0:ip1])
      data_on_grid[ip0:ip1] = dot(ao, data)
  elif 'den' in dtype:
    for ip0, ip1 in lib.prange(0, ngrids, blksize):
      if verbose: print(f"Done with {ip0}/{ngrids} = {ip0/ngrids:0.2f}.")
      ao = mol.eval_gto(GTOval, grid[ip0:ip1])
      data_on_grid[ip0:ip1] = eval_rho(mol, ao, data)
  else:
    raise AssertionError("dtype={dtype} not understood, try 'orbital' or density'.")

  return data_on_grid

def make_grid(voxel=eye(3), npoints=(10,10,10), origin=(0,0,0), ptshift=(0,0,0)):
  ''' Make a grid of points defining a set of voxels.
  Args:
    voxel (array): defines the spaces between points.
    npoints (tuple): the number of voxels repeated in each dimension.
    origin (array-like): the corner of the volume defined by the first voxel. 
    ptshift (tuple): Shift grid from the origin by this many voxels. For example, use npoints/2 for centering the origin.
  Returns:
    grid (array): Points of the grid in cartesian coordinates. 
    voxgrid (array): Points of the grid in the basis of the voxel vectors.
  '''
 
  voxgrid = asarray(meshgrid(
    range(-ptshift[0],npoints[0]-ptshift[0]), 
    range(-ptshift[1],npoints[1]-ptshift[1]), 
    range(-ptshift[2],npoints[2]-ptshift[2]), 
    indexing='ij'))
  voxgrid = voxgrid.reshape(voxgrid.shape[0], voxgrid.size//voxgrid.shape[0]).T
  grid = voxgrid @ voxel + asarray(origin)

  return grid, voxgrid

# I think this is obsolete with make_grid.
def make_grid_slow(voxel=eye(3),npoints=(10,10,10),origin=(0,0,0),skip=((0,0),(0,0),(0,0))):
  ''' Make a grid at npoints multiples of voxvecs.'''
  size = ((npoints[0]-sum(skip[0])),(npoints[1]-sum(skip[1])),(npoints[2]-sum(skip[2])))
  grid = empty((size[0],size[1],size[2],3))

  for i in range(size[0]):
    ipart = (i+skip[0][0]+origin[0]*npoints[0])%npoints[0] * voxel[0]
    for j in range(size[1]):
      ijpart = (j+skip[1][0]+origin[1]*npoints[1])%npoints[1] * voxel[1] + ipart
      for k in range(size[2]):
        grid[i,j,k,:] = (k+skip[2][0]+origin[2]*npoints[2])%npoints[2] * voxel[2] + ijpart

  return grid.reshape(size[0]*size[1]*size[2],3)

def _set_nonzero_origin(cart, latvecs, origin):
  frac = inv(latvecs) @ cart
  frac = (frac + origin[:,None])%1.0

  return latvecs @ frac 

if __name__=='__main__':
  main()
