#!/usr/bin/python
from numpy  import array,dot
from mython import lines2str,gen_qsub
from os import getcwd

# Reads a CRYSTAL input file. Can easily be modified to get more information.
def read_cryinp(inpf):
  inpstr = ''
  res = {}
  titleline = True
  for line in inpf:
    if titleline: titleline=False; continue
    inpstr += line

  lines = inpstr.split('\n')
  pos = 2

  # Read in geometry section which is primarily position-based.
  res['group']    = int(lines[pos]);                pos += 1
  res['latparms'] = map(float,lines[pos].split());  pos += 1
  res['natoms']   = int(lines[pos]);                pos += 1 #TODO Not always correct!
  res['atypes'], res['apos'] = [],[]
  for i in range(res['natoms']):
    line = lines[pos].split()
    res['atypes'].append(int(line[0]))
    res['apos'].append(map(float,line[1:]))
    pos += 1

  # Rest of input not position based, and may not be separated correctly by
  # newlines. TODO: this might be true of geometry as well.

  inpl = inpstr.split()

  # Now read in important keyword arguements.
  res['calculation'] = 'hf'
  res['effchg'] = []
  res['supercell'] = [1,0,0,0,1,0,0,0,1]
  res['tolinteg'] = [6,6,6,6,12]
  while pos < len(inpl):
    if 'SUPERCELL' == inpl[pos]:
      res['supercell'] =  map(int,inpl[pos+1:pos+10])
      pos += 11
      continue

    if 'DFT' == inpl[pos]:
      res['calculation'] = 'dft'
      pos += 1
      continue

    if 'CORRELAT' == inpl[pos]:
      res['correlation'] = inpl[pos+1].lower()
      pos += 2
      continue

    if 'EXCHANGE' == inpl[pos]:
      res['exchange'] = inpl[pos+1].lower()
      pos += 2
      continue

    if 'HYBRID' == inpl[pos]:
      res['mixing'] = float(inpl[pos+1])
      pos += 2
      continue

    if 'PBE0' == inpl[pos]:
      res['correlation'] = 'pbe'
      res['exchange'] = 'pbe'
      res['mixing'] = 25.0
      pos += 1
      continue

    # This currently depends on the order of the entry when it doesn't have to.
    if 'INPUT' == inpl[pos]:
      res['effchg'].append(float(inpl[pos+1]))
      pos += 1
      continue

    if 'SHRINK' == inpl[pos]:
      res['kdens'] = int(inpl[pos+1].split()[0])
      pos += 2
      continue

    if 'TOLINTEG' == inpl[pos]:
      res['tolinteg'] = int(inpl[pos+1].split()[0])
      pos += 2
      continue

    if 'TOLDEE' == inpl[pos]:
      res['tole'] = int(inpl[pos+1])
      pos += 2
      continue

    if 'SPINLOCK' == inpl[pos]:
      res['spinlock'] = int(inpl[pos+1].split()[0])
      pos += 2
      continue

    if 'FMIXING' == inpl[pos]:
      res['fmixing'] = int(inpl[pos+1])
      pos += 2
      continue

    if 'BROYDEN' == inpl[pos]:
      res['broyden'] = [float(inpl[pos+1])] + map(int,inpl[pos+2:pos+4])
      pos += 5
      continue

    if 'ATOMSPIN' == inpl[pos]:
      nspin = int(inpl[pos+1])
      spinl = inpl[pos+2:pos+2+2*nspin]
      res['atomspin'] = [(int(spinl[i]),int(spinl[i+1])) for i in range(nspin)]

    pos += 1
  return res

def read_cryout(inpf):
  res = {}
  line = 'start'
  spins = []
  while line != '':
    line = inpf.readline()
    if 'TOTAL ATOMIC SPINS' in line:
      spins = []
      line  = inpf.readline()
      while ('TTT' not in line) and (line != ''):
        spins += map(float,line.split())
        line   = inpf.readline()

    if 'SCF ENDED' in line:
      spl = line.split()
      if spl[4] != 'CONVERGENCE':
        print "Severe warning: DFT SCF not converged!"
      res['dft_energy'] = float(spl[8])
  if spins != []:
    res['dft_moments'] = spins
  return res

def gen_properties(cryinp,natoms,kpath,denom,projs,
                   get_ef=1,nprint=0,title='blank',
                   above=5,below=8,
                   kresfactor=20,eres=200,
                   npoly=25,):
  """ User-friendly properties input generator, for bandstructure and DOS.

      Factor "denom" out of kpath to make it integer.
      TODO: negative DOSS projections (onto atoms)
  """
  inp = read_cryinp(cryinp)
  shrink = inp['kdens']
  epera  = array(inp['effchg'])
  
  # This may only make sense when the spin channels are equal.
  filled = int(round(dot(array(natoms),epera))) / 2

  newklines = [['NEWK']]
  newklines.append([shrink,shrink*2])
  newklines.append([get_ef,nprint])

  nline       = len(kpath)
  kres        = kresfactor*nline
  highest     = filled + above
  lowest      = filled - below
  outtofile   = 1
  printevals  = 0
  bandlines = [['BAND']]
  bandlines.append([title])
  bandlines.append([nline,denom,kres,lowest,highest,outtofile,printevals])
  for si in range(len(kpath)-1):
    bandlines.append(kpath[si]+[' ']+kpath[si+1])
  bandlines.append(kpath[-1]+[' ']+kpath[0])

  nproj      =  len(projs)
  printopts  = 0
  doslines = [['DOSS']]
  doslines.append([nproj,eres,lowest,highest,outtofile,npoly,printopts])
  for pi in range(len(projs)):
    doslines.append([len(projs[pi]),' ']+projs[pi])

  outlines = newklines + bandlines + newklines + doslines + [["END"]]
  return lines2str(outlines)

def gen_cry2qwalk(dftfn):
  loc = getcwd()
  root = dftfn.replace('.d12','')
  propfn = root+'.d3'
  dftdat = read_cryinp(open(dftfn,'r'))
  shrink = dftdat['kdens']
  proplines = []
  proplines.append('NEWK')
  proplines.append('{0} {1}'.format(shrink,2*shrink))
  proplines.append('0 1')
  proplines.append('67 999')
  proplines.append('END')
  with open(propfn,'w') as propf:
    propf.write('\n'.join(proplines))
  exe = '~/bin/properties < {0}'.format(propfn)
  out = propfn+'.out'
  fc = []
  fc.append('cat {root}.d12.out {root}.d3.out > tmp'.format(root=root))
  fc.append('addline=`grep -n "HAMILTONIAN EIGENVECTORS" tmp | head -1 | cut -d":" -f1`')
  fc.append('sed "$addline i NEWK EIGENVECTORS" tmp > crystal2qmc.inp')
  fc.append('rm -f tmp')
  fc.append('~/bin/crystal2qmc -c -o {root} crystal2qmc.inp &>> crystal2qmc.out'.format(root=root))
  gen_qsub(exe,stdout=out,loc=loc,
           name='cvt_'+root,
           time='02:00:00',
           queue='secondary',
           final_commands=fc)

