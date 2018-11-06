import sys
import numpy as np
from numpy.linalg import det
from scipy.interpolate import griddata
from scipy.signal import butter, lfilter
from copy import deepcopy

#####################################
def read_cube(inpf,qwalk_patch=False):
  if type(inpf)==str: inpf = open(inpf,'r')
  cube={}
  cube['comment']=inpf.readline()
  cube['type']=inpf.readline()
  spl=inpf.readline().split()
  #cube['natoms']=int(spl[0])
  cube['natoms']=round(float(spl[0]))
  if qwalk_patch:
    cube['origin']=np.zeros(3)
  else:
    cube['origin']=np.array(spl[1:],dtype=float)
  cube['ints']=np.array([0,0,0])
  cube['latvec']=np.zeros((3,3))
  for i in range(0,3):
    spl=inpf.readline().split()
    cube['ints'][i]=int(spl[0])
    cube['latvec'][i,:]=np.array(spl[1:],dtype=float)
  natoms=cube['natoms']
  cube['atomname']=[]
  cube['atomxyz']=np.zeros((natoms,3))
  for i in range(0,natoms):
    spl=inpf.readline().split()
    cube['atomname'].append(spl[0])
    cube['atomxyz'][i,:]=np.array(spl[2:],dtype=float)
  cube['data']=np.zeros(cube['ints'])
  vector=[]
  while True:
    spl=inpf.readline().split()
    if len(spl) < 1:
      break
    vector.extend(map(float,spl))
  nread=len(vector)
  count=0
  for x in range(0,cube['ints'][0]):
    for y in range(0,cube['ints'][1]):
      for z in range(0,cube['ints'][2]):
        cube['data'][x,y,z]=vector[count]
        count+=1
        #if count >= nread:
        #  break;
      #if count>=nread:
      #  break
    #if count >= nread:
    #  break

  return cube

#####################################
def write_cube(cube, outf):
  if type(outf)==str: outf = open(outf,'w')
  outf.write(cube['comment'])
  outf.write(cube['type'])
  outf.write(str(cube['natoms'])+" {} {} {}".format(*cube['origin']))
  outf.write("\n")
  for i in range(0,3):
    outf.write("%i "%cube['ints'][i])
    outf.write(" %g %g %g \n"%(cube['latvec'][i,0],cube['latvec'][i,1],cube['latvec'][i,2]))
  natoms=cube['natoms']
  for i in range(0,natoms):
    outf.write("%s 0.0 "%cube['atomname'][i])
    outf.write(" %g %g %g \n"%(cube['atomxyz'][i,0],cube['atomxyz'][i,1],cube['atomxyz'][i,2]))
  count=0
  for x in range(0,cube['ints'][0]):
    for y in range(0,cube['ints'][1]):
      for z in range(0,cube['ints'][2]):
        outf.write("%g "%cube['data'][x,y,z])
        count+=1
        if count%5==0:
          outf.write('\n')
  outf.write('\n')

#####################################
def write_xsf(cube,outf):
  if type(outf)==str: outf=open(outf,'w')
  outf.write("CRYSTAL\n")
  outf.write("PRIMVEC\n")
  natoms=cube['natoms']
  for i in range(0,3):
    npts=cube['ints'][i]
    outf.write(" %g %g %g \n"%(npts*cube['latvec'][i,0],npts*cube['latvec'][i,1],npts*cube['latvec'][i,2]))
  outf.write("PRIMCOORD\n")
  outf.write("%i 1\n"%natoms)
  for i in range(0,natoms):
    outf.write("%s "%cube['atomname'][i])
    outf.write(" %g %g %g \n"%(cube['atomxyz'][i,0],cube['atomxyz'][i,1],cube['atomxyz'][i,2]))
  outf.write("BEGIN_BLOCK_DATAGRID_3D\n cube_file_conversion \n")
  outf.write("BEGIN_DATAGRID_3D\n")
  outf.write("%i %i %i\n"%(cube['ints'][0],cube['ints'][1],cube['ints'][2]))
  outf.write("0.0 0.0 0.0\n")
  for i in range(0,3):
    npts=cube['ints'][i]
    outf.write(" %g %g %g \n"%(npts*cube['latvec'][i,0],npts*cube['latvec'][i,1],npts*cube['latvec'][i,2]))
  
  count=0
  for z in range(0,cube['ints'][2]):
    for y in range(0,cube['ints'][1]):
      for x in range(0,cube['ints'][0]):
        outf.write("%g "%cube['data'][x,y,z])
        count+=1
        if count%5==0:
          outf.write('\n')
  outf.write('\n')
  outf.write("END_DATAGRID_3D\n")
  outf.write("END_BLOCK_DATAGRID_3D\n")
  
#####################################
def integrate(cube):
  """Numerically integrate the density.
  
  Appoximates integral by simple sum."""
  vol=abs(det(cube['latvec']))
  return np.sum(cube['data'])*vol

#####################################
def integrate_abs(cube):
  """Numerically integrate the absolute value of the density.
  
  Appoximates integral by simple sum."""
  vol=abs(det(cube['latvec']))
  return np.sum(abs(cube['data']))*vol

#####################################
def cabs(cube):
  """Take absolute value of cube data."""
  newcube = deepcopy(cube)
  newcube['data']=abs(newcube['data'])
  return newcube

#####################################
def normalize_abs(cube,Nelec=1):
  """Normalize the density so the integral over all space yeilds Nelec.
  
  Appoximates integral by simple sum."""
  vol=abs(det(cube['latvec']))
  norm=np.sum(abs(cube['data']))*vol
  cube['data']*=(float(Nelec)/norm)
  return cube

#####################################
def freq_cutoff(cube,freq_cutoff=0.90):
  """ Cutoff frequencies of the signal with size freq_cutoff * the
  maximum or higher, in place."""
  # Frequency distribution is fastest in middle, slowest at end, and positive in
  # first half.
  max_val = cube['data'].max()

  fft = np.fft.fftn(cube['data'])
  endfreqs = np.array([int(round(s*freq_cutoff/2.)) for s in fft.shape])

  print("Cutting off ",)
  for si,s in enumerate(fft.shape): print("%d "%(s - 2*int(round(endfreqs[si]))),)
  print("frequencies, out of ",)
  print("%d %d %d"%fft.shape)

  for dim in range(fft.ndim):
    fft = fft.swapaxes(dim,fft.ndim-1)
    for d1 in fft:
      for d2 in d1:
        d2[endfreqs[dim]:-endfreqs[dim]+1] = 0.2*d2[endfreqs[dim]:-endfreqs[dim]+1]
    fft = fft.swapaxes(dim,fft.ndim-1)
  cube['data'] = np.fft.ifftn(fft)
  if abs(cube['data'].imag).max() > 1e-16:
    print("Warning, inverting FFT may not be completely real!")
  cube['data'] = cube['data'].real
  cube['data'] *= (max_val / cube['data'].max())

#####################################
def butter_cutoff(cube,crit_freq=1,order=4):
  """
  Simple wrapper for scipy routines to perform Butterworth low pass filter.
  
  crit_freq = point where gain drops to 1/sqrt(2) of passband. 1 is defined as
  the Nyquist frequency.
  order = Order of Butterworth function, which controls steepness.
  """
  b,a = butter(order, crit_freq)
  cube['data'] = lfilter(b,a,cube['data'])
  #cube['data'] = lfilter(b,a,cube['data'],0)
  #cube['data'] = lfilter(b,a,cube['data'],1)
  #cube['data'] = lfilter(b,a,cube['data'],2)

#####################################
def gaussian_averager(cube,sigma=3,nbr_dist=1,repeat=1):
  """ Average each point in the cube file with blob_range neighbors in each
  direction, weighted by a Gaussian with SD sigma."""

  nd = nbr_dist
  total_steps = 0.
  for ii in range(-nd,nd+1):
    for jj in range(-(nd-abs(ii)),(nd-abs(ii))+1):
      for kk in range(-(nd-abs(ii)-abs(jj)),(nd-abs(ii)-abs(jj))+1):
        total_steps += 1
  total_steps *= repeat
  done_steps = 0.

  for iteration in range(repeat):
    new = np.zeros(cube['data'].shape)
    def wf(x):
      return np.exp(-x**2/(2*sigma**2))

    for ii in range(-nd,nd+1):
      wi = wf(ii)
      for jj in range(-(nd-abs(ii)),(nd-abs(ii))+1):
        wj = wf(jj)
        for kk in range(-(nd-abs(ii)-abs(jj)),(nd-abs(ii)-abs(jj))+1):
          wk = wf(kk)
          print("Finished {:5.2%}. ".format(done_steps/total_steps))
          for i in range(cube['data'].shape[0]):
            ip = (i+ii)%cube['data'].shape[0]
            for j in range(cube['data'].shape[1]):
              jp = (j+jj)%cube['data'].shape[1]
              for k in range(cube['data'].shape[2]):
                kp = (k+kk)%cube['data'].shape[2]
                new[i,j,k] += cube['data'][ip,jp,kp]*wi*wj*wk
          done_steps += 1
    cube['data'] = new

#####################################
def sub_cubes(poscube,negcube):
  """Subtract two cube files.
  Note: you may need to normalize these appropriately first."""

  subcube = deepcopy(poscube)
  subcube['data'] -= negcube['data']
  return subcube

#####################################
def add_cubes(cube1,cube2,N1=1,N2=1):
  """Add two cube files.
  Note: you may need to normalize these appropriately first."""
  addcube = deepcopy(cube1)
  addcube['data'] += cube2['data']
  #addcube['data'] /= abs(addcube['data']).sum()
  return addcube

#####################################
def mul_cubes(cube1,cube2,N1=1,N2=1):
  """Multiply two cube files, pointwise.
  Note: you may need to normalize these appropriately first."""
  mulcube = deepcopy(cube1)
  mulcube['data'] *= cube2['data']
  #mulcube['data'] /= abs(mulcube['data']).sum()
  return mulcube

#####################################
# Used for interpolation scheme
def nearest(point,cube):
  """Find the value in the cube file located closest to point."""
  a = np.array([ np.dot(cube['latvec'][i],point)/np.dot(cube['latvec'][i],cube['latvec'][i])
              for i in range(3) ]).round()
  #print a % cube['ints']
  return cube['data'][tuple(map(int,a % cube['ints']))]

#####################################
# Used for interpolation scheme
def linear(point,cube):
  """Compute the linear extrapolation to the point using closest available
  points."""
  latvec = cube['latvec'] 
  ints = cube['ints']
  # point in the basis of latvec.
  pnb = np.array([ np.dot(latvec[i],point)/np.dot(latvec[i],latvec[i])
              for i in range(3) ])
  # Round up and down to get points on the lattice.
  neighbors = [(ax,ay,az)
    for ax in map(int,[np.floor(pnb[0]),np.ceil(pnb[0])])
    for ay in map(int,[np.floor(pnb[1]),np.ceil(pnb[1])]) 
    for az in map(int,[np.floor(pnb[2]),np.ceil(pnb[2])])]
  vals = np.array([cube['data'][tuple(n)] for n in (neighbors%ints)])
  vals = vals.reshape(2,2,2)
  tvol = abs(det(latvec))
  vols = np.array([abs(det((n-pnb)*latvec)) for n in neighbors]).reshape(2,2,2)
  wght = np.zeros(vals.shape)
  # Weights for average are the subvolumes of opposing corner.
  for i in range(2):
    for j in range(2):
      for k in range(2):
        wght[i,j,k] = vols[1-i,1-j,1-k]/tvol
  return np.sum(wght*vals)

#####################################
def interp_cube(cube, pos, res=(10,10), method='nearest', atrad=0.0):
  """Interpolate cube in plane defined by three points, pos, with res points, using
  method to interpolate, and ensuring atrad radius around each atom is included.
  By design, tries to put the longer axis on the x axis. To change this, you'd
  need to switch pvm and pvo."""
  latvec = cube['latvec']
  ints = cube['ints']
  
  # Classify points to show all atoms while minimizing extra space.
  pidx = (np.array((pos[1]-pos[2], pos[2]-pos[0], pos[0]-pos[1]))**2).sum(axis=1).argsort()
  orig = pos[pidx[ 1]]
  pvm  = pos[pidx[-1]] - orig # Main plotting axis (static).
  pvo  = pos[pidx[ 0]] - orig # Orthogonal plotting axis.
  # Idea is that normally you'd want to orthogonalize the sortest axis.
  # Moving the axis is what creates extra space in the plot.

  # Orthogonalize orth. axis, since this is how most plots are.
  pvo -= np.sum(pvo*pvm)/np.sum(pvm**2) * pvm
  # Make it a square plot.
  #pvo *= (sum(pvm**2)/sum(pvo**2))**.5

  # Add buffer to contain all of atom. Extra space in each plot dir.
  buffm = atrad/(np.sum(pvm**2))**.5 * pvm 
  buffo = atrad/(np.sum(pvo**2))**.5 * pvo

  pvm += 2*buffm
  pvo += 2*buffo

  e1 = pvm / float(res[0]-1)
  e2 = pvo / float(res[1]-1)

  # Domain of output.
  odom = np.array( [ i*e1 + j*e2
                   for i in range(res[0])
                   for j in range(res[1]) ])
  odom += orig - buffm - buffo
  #odom.shape = (cumprod(odom.shape[:2])[-1], odom.shape[2])

  # Wrap at zone boundaries.
  basis = latvec * ints
  for j in range(len(odom)):
    odom[j] = np.sum( [(np.dot(basis[i],odom[j])/np.dot(basis[i],basis[i]))%1*basis[i]
                    for i in range(3)], axis=1 )

  # Compute the closest point of the atoms to the plane, and how far away they
  # are.
  atpos = [a - orig + buffm + buffo for a in cube['atomxyz']]

  # includes periodic images of neighboring cells.
  atpos += [a-b for b in basis for a in atpos] + [a+b for b in basis for a in atpos]


  acoor = [( np.dot(pvm,a) / np.dot(pvm,pvm)**.5,
             np.dot(pvo,a) / np.dot(pvo,pvo)**.5 )
           for a in atpos]
  adist = [np.dot( np.cross(pvm,pvo)/(np.dot(pvm,pvm)*np.dot(pvo,pvo))**.5, a ) 
           for a in atpos]

  X = np.linspace(0,np.sum(pvm**2)**.5,res[0])
  Y = np.linspace(0,np.sum(pvo**2)**.5,res[1])

  if method=='nearest':
    Z = np.array([nearest(point,cube) for point in odom]).reshape(res)
  elif method=='linear':
    Z = np.array([linear(point,cube) for point in odom]).reshape(res)
  else:
    print('Interpolation method is not implemented yet.')
  
  return {'points':(X,Y), 'data':Z, 'acoor':np.array(acoor), 'adist':np.array(adist)}

if __name__=="__main__":
  import sys

  implemented = ['add','sub','xsf','integrate_abs']

  if len(sys.argv) < 2:
    raise AssertionError("""
    Usage: python cube.py <operation> <list of input cubes> 
    Operations: {}.""".format(implemented))

  if sys.argv[1] == "add":
    if len(sys.argv) != 4:
      raise AssertionError("Add needs two input cubes.")
    outfn = input("Output cube: ")
    with open(outfn,'w') as outf:
      write_cube(
          add_cubes(
              read_cube(open(sys.argv[2],'r')),
              read_cube(open(sys.argv[3],'r'))
            ),
          outf
        )
  elif sys.argv[1] == "sub":
    if len(sys.argv) != 4:
      raise AssertionError("Subtract needs two input cubes.")
    outfn = input("Output cube: ")
    with open(outfn,'w') as outf:
      write_cube(
          sub_cubes(
              read_cube(open(sys.argv[2],'r')),
              read_cube(open(sys.argv[3],'r'))
            ),
          outf
        )
  elif sys.argv[1] == "xsf":
    if len(sys.argv) != 3:
      raise AssertionError("Converter needs exactly one input cube.")
    outfn = input("Output xsf: ")
    with open(outfn,'w') as outf:
      write_xsf(
          read_cube(open(sys.argv[2],'r')),
          outf
        )
  elif sys.argv[1] == "integrate_abs":
    if len(sys.argv) != 3:
      raise AssertionError("Integration needs exactly one input cube.")
    print("Abs. integral: ",integrate_abs(read_cube(open(sys.argv[2],'r'))))

  else:
    raise AssertionError("""
        Sorry, '{}' keyword isn't implemented.
        Implemented operations are {}.
        It's probably trivial to add this operation yourself, you should do it and
        push the result!""".format(sys.argv[1],implemented)
      )
