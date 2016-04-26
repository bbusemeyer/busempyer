import sys
from numpy import array,linspace,zeros,ones,cumprod,floor,ceil,dot,cross,sum,nan_to_num,errstate
import numpy as np
from numpy.linalg import det
from scipy.interpolate import griddata
from scipy.signal import butter, lfilter
from copy import deepcopy

#####################################
def read_cube(inpf):
  cube={}
  cube['comment']=inpf.readline()
  cube['type']=inpf.readline()
  spl=inpf.readline().split()
  cube['natoms']=int(spl[0])
  cube['origin']=map(float, spl[1:])
  cube['ints']=array([0,0,0])
  cube['latvec']=zeros((3,3))
  for i in range(0,3):
    spl=inpf.readline().split()
    cube['ints'][i]=int(spl[0])
    cube['latvec'][i,:]=map(float,spl[1:])
  natoms=cube['natoms']
  cube['atomname']=[]
  cube['atomxyz']=zeros((natoms,3))
  for i in range(0,natoms):
    spl=inpf.readline().split()
    cube['atomname'].append(spl[0])
    cube['atomxyz'][i,:]=map(float,spl[2:])
  cube['data']=zeros(cube['ints'])
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
  outf.write(cube['comment'])
  outf.write(cube['type'])
  outf.write(' '.join(map(str,[cube['natoms']]+cube['origin'])))
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
  return sum(cube['data'])*vol

#####################################
def integrate_abs(cube):
  """Numerically integrate the absolute value of the density.
  
  Appoximates integral by simple sum."""
  vol=abs(det(cube['latvec']))
  return sum(abs(cube['data']))*vol

#####################################
def normalize_abs(cube,Nelec=1):
  """Normalize the density so the integral over all space yeilds Nelec.
  
  Appoximates integral by simple sum."""
  vol=abs(det(cube['latvec']))
  norm=sum(abs(cube['data']))*vol
  cube['data']*=(float(Nelec)/norm)

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
def sub_cubes(poscube,negcube,Npos=1,Nneg=1):
  """Subtract two cube files.

  Will normalize first, which requres Nup and Ndn arguments if the number of
  electrons in the two cubes differs."""
  normalize_abs(poscube,Npos)
  normalize_abs(negcube,Nneg)
  subcube = deepcopy(poscube)
  subcube['data'] -= negcube['data']
  return subcube

#####################################
def add_cubes(cube1,cube2,N1=1,N2=1):
  """Add two cube files.

  Will normalize first, which requres N1 and N2 arguments if the number of
  electrons in the two cubes differs."""
  normalize_abs(cube1,N1)
  normalize_abs(cube2,N2)
  addcube = deepcopy(cube1)
  addcube['data'] += cube2['data']
  addcube['data'] /= abs(addcube['data']).sum()
  return addcube

#####################################
# Used for interpolation scheme
def nearest(point,cube):
  """Find the value in the cube file located closest to point."""
  a = array([ dot(cube['latvec'][i],point)/dot(cube['latvec'][i],cube['latvec'][i])
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
  pnb = array([ dot(latvec[i],point)/dot(latvec[i],latvec[i])
              for i in range(3) ])
  # Round up and down to get points on the lattice.
  neighbors = [(ax,ay,az)
    for ax in map(int,[floor(pnb[0]),ceil(pnb[0])])
    for ay in map(int,[floor(pnb[1]),ceil(pnb[1])]) 
    for az in map(int,[floor(pnb[2]),ceil(pnb[2])])]
  vals = array([cube['data'][tuple(n)] for n in (neighbors%ints)])
  vals = vals.reshape(2,2,2)
  tvol = abs(det(latvec))
  vols = array([abs(det((n-pnb)*latvec)) for n in neighbors]).reshape(2,2,2)
  wght = zeros(vals.shape)
  # Weights for average are the subvolumes of opposing corner.
  for i in range(2):
    for j in range(2):
      for k in range(2):
        wght[i,j,k] = vols[1-i,1-j,1-k]/tvol
  return sum(wght*vals)

#####################################
def interp_cube(cube, pos, res=(10,10), method='nearest', atrad=0.0):
  """Interpolate cube in plane defined by three points, pos, with res points, using
  method to interpolate, and ensuring atrad radius around each atom is included.
  By design, tries to put the longer axis on the x axis. To change this, you'd
  need to switch pvm and pvo."""
  latvec = cube['latvec']
  ints = cube['ints']
  
  # Classify points to show all atoms while minimizing extra space.
  pidx = (array((pos[1]-pos[2], pos[2]-pos[0], pos[0]-pos[1]))**2).sum(axis=1).argsort()
  orig = pos[pidx[ 1]]
  pvm  = pos[pidx[ 0]] - orig # Main plotting axis (static).
  pvo  = pos[pidx[-1]] - orig # Orthogonal plotting axis.

  # Orthogonalize orth. axis, since this is how most plots are.
  pvo -= sum(pvo*pvm)/sum(pvm**2) * pvm
  # Make it a square plot.
  #pvo *= (sum(pvm**2)/sum(pvo**2))**.5

  # Add buffer to contain all of atom. Extra space in each plot dir.
  buffm = atrad/(sum(pvm**2))**.5 * pvm 
  buffo = atrad/(sum(pvo**2))**.5 * pvo

  pvm += 2*buffm
  pvo += 2*buffo

  e1 = pvm / float(res[0]-1)
  e2 = pvo / float(res[1]-1)

  # Domain of output.
  odom = array( [ i*e1 + j*e2
                   for i in range(res[0])
                   for j in range(res[1]) ])
  odom += orig - buffm - buffo
  #odom.shape = (cumprod(odom.shape[:2])[-1], odom.shape[2])

  # Wrap at zone boundaries.
  basis = latvec * ints
  for j in range(len(odom)):
    odom[j] = sum( [(dot(basis[i],odom[j])/dot(basis[i],basis[i]))%1*basis[i]
                    for i in range(3)], axis=1 )

  # Compute the closest point of the atoms to the plane, and how far away they
  # are.
  atpos = [a - orig + buffm + buffo for a in cube['atomxyz']]

  # includes periodic images of neighboring cells.
  atpos += [a-b for b in basis for a in atpos] + [a+b for b in basis for a in atpos]


  acoor = [( dot(pvm,a) / dot(pvm,pvm)**.5,
             dot(pvo,a) / dot(pvo,pvo)**.5 )
           for a in atpos]
  adist = [dot( cross(pvm,pvo)/(dot(pvm,pvm)*dot(pvo,pvo))**.5, a ) 
           for a in atpos]

  X = linspace(0,sum(pvm**2)**.5,res[0])
  Y = linspace(0,sum(pvo**2)**.5,res[1])

  if method=='nearest':
    Z = array([nearest(point,cube) for point in odom]).reshape(res)
  elif method=='linear':
    Z = array([linear(point,cube) for point in odom]).reshape(res)
  else:
    print('Interpolation method is not implemented yet.')
  
  return {'points':(X,Y), 'data':Z, 'acoor':array(acoor), 'adist':array(adist)}
