import sys
from numpy import array,linspace,zeros,ones,cumprod,floor,ceil,dot,cross,sum,nan_to_num,errstate
from numpy.linalg import det
from scipy.interpolate import griddata

def read_cube(inpf):
  cube={}
  cube['comment']=inpf.readline()
  cube['type']=inpf.readline()
  spl=inpf.readline().split()
  #cube['natoms']=int(spl[0])
  cube['natoms']=map(int,spl)
  cube['origin']=map(float, spl[1:])
  cube['ints']=array([0,0,0])
  cube['latvec']=zeros((3,3))
  for i in range(0,3):
    spl=inpf.readline().split()
    cube['ints'][i]=int(spl[0])
    cube['latvec'][i,:]=map(float,spl[1:])
  natoms=cube['natoms'][0]
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

def write_cube(cube, outf):
  outf.write(cube['comment'])
  outf.write(cube['type'])
  for i in cube['natoms']:
    outf.write(" %i "%i)
  f.write("\n")
  for i in range(0,3):
    outf.write("%i "%cube['ints'][i])
    outf.write(" %g %g %g \n"%(cube['latvec'][i,0],cube['latvec'][i,1],cube['latvec'][i,2]))
  natoms=cube['natoms'][0]
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
  
def normalize_abs(cube):
  vol=abs(linalg.det(cube['latvec']))
  norm=sum(abs(cube['data']))*vol
  cube['data']/=norm
  return cube

# Used for interpolation scheme
def nearest(point,cube):
  a = array([ dot(cube['latvec'][i],point)/dot(cube['latvec'][i],cube['latvec'][i])
              for i in range(3) ]).round()
  #print a % cube['ints']
  return cube['data'][tuple(map(int,a % cube['ints']))]

def linear(point,cube):
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

# Interpolate cube in plane defined by three points, pos, with res points, using
# method to interpolate, and ensuring atrad radius around each atom is included.
# By design, tries to put the longer axis on the x axis. To change this, you'd
# need to switch pvm and pvo.
def interp_cube(cube, pos, res=(10,10), method='nearest', atrad=0.0):
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
    print 'Interpolation method is not implemented yet.'
  
  return {'points':(X,Y), 'data':Z, 'acoor':array(acoor), 'adist':array(adist)}
