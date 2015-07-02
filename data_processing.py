#!/usr/bin/python
from numpy            import array,linalg
from mython           import Ldict
from subprocess       import call
from dm_tools         import read_dm
from read_numberfluct import read_number_dens,moments
from cryfiles_io      import read_cryinp, read_cryout
from os               import getcwd
import sys

# Reads a qwalk input section.
def read_section(inp,key,pos):
  res = Ldict()
  while inp[pos] != '}':
    if isinstance(inp[pos],float):           # If it's a number.
      while isinstance(inp[pos],float):
        #print 'Found data',key,inp[pos]
        res.append(key,inp[pos])
        pos += 1
    elif inp[pos] == '{':                    # Else if it's demarking a section,
      if isinstance(inp[pos+1],str):         # which could have keywords,
        label = inp[pos+1].lower()
        pos += 2
        #print 'Reading section',key,label
        val,pos = read_section(inp,key,pos)
        if label != False:
          val['label'] = label
        res.append(key,val)
      else:                                  # or just numbers.
        pos += 1
        val = []
        while isinstance(inp[pos],float):
          val.append(inp[pos])
          pos += 1
        if len(val) == 1: val = val[0]
        #print 'Found data',key,val
        res[key] = val
        pos += 1
    else:                                    # Else it's a keyword.
      key = inp[pos].lower()
      if key not in res.keys():
        #print 'Setting',key
        res[key] = True
      pos += 1
  pos += 1
  return res,pos

# Reads a qwalk input file.
def read_qfile(inpf):
  inpstr = ''
  for line in inpf:
    if '#' in line: # TODO: Needs to be fixed when '#' isn't the first thing in line.
      print 'Warning, reading commented lines is incomplete!'
      continue
    inpstr += line
  # Ensure correct splitting. This is inefficient for large files.
  inpstr = inpstr.replace('{',' { ')
  inpstr = inpstr.replace('}',' } ')
  inp    = inpstr.split() + ['}'] # Now I can use "read_section" on the file!
  for i in range(len(inp)): # Make everything numeric into floats.
    try:                inp[i] = float(inp[i])
    except ValueError:  pass
  return read_section(inp,inpf.name,0)[0]

# Temporary function to convert file names to metadata about the calculations.
# In the future, an optional metadata file should contain overall qualitiative
# descriptions of the data, like magnetic ordering, or pressure. These things
# are redundant, but convenient for interpreting sets of numbers together.
def convert_to_metadata(froot):
  mconv = {'che':'checkerboard',
           'str':'collinear',
           'bic':'bicollinear',
           'fst':'collinear, flip 1',
           'dim':'collinar, flip 2',
           'sta':'staggered'}
  basename = froot.split('/')[-1]
  mag = mconv[basename[:3]]
  prs = float(basename[3:].translate(None,'_prs'))
  return {'mag':mag,'prs':prs}

# Use golsing to read out average energy and error.
def read_qenergy(logfile,gosling='./gosling'):
  statfilen = logfile.name.replace('.log','.stat')
  with open(statfilen,'w') as out:
    try:
      call([gosling, logfile.name], stdout = out)
    except OSError:
      print "Cannot find gosling!"
      exit()
  with open(statfilen,'r') as F:
    for line in F:
      if 'total_energy0' in line:
        spl = line.split()
        return {'egy':spl[1],'err':spl[3],'var':spl[5]}
  print 'ERROR: cannot find total_energy0 in stat file.'
  return {'egy':None,'err':None,'var':None}

# Read a directory which has multiple files with data into a single dictionary
# with relevant information.
def read_dir_forlucas(froot,gosling='./gosling'):
  """ Reads a CRYSTAL + QWalk directory's data into a dictionary 
  
  forlucas implies its the keys that he specified via email."""

  ############################################################################
  # This first section should be edited to reflect naming conventions!
  dftfile   = froot+'.d12'
  # Start by taking only the real k-points, since I'm sure these are on solid
  # ground, and have enough sample points.
  realk = array([1,3,8,10,27,29,34,36]) - 1
  ############################################################################

  bres = {} # Data that is common to all k-points.
  ress = {} # Dict of all k-point data in directory.

  metad = convert_to_metadata(froot)
  bres['authors']  = 'BL'
  bres['ordering'] = metad['mag']

  try:
    inpf = open(dftfile,'r')
    dftdat = read_cryinp(inpf)
    bres['a'] = dftdat['latparms'][0]
    bres['c'] = dftdat['latparms'][1]
    bres['se_height'] = dftdat['apos'][dftdat['atypes'].index(234)][-1]
    bres['hybrid_mixing'] = dftdat['mixing']
  except IOError:
    print "There's no dft in this directory!"
    return {}

  for rk in realk:
    kroot = froot + '_' + str(rk)

    print "Working on",kroot+"..." 
    try:
      inpf = open(kroot+'.sys','r')
      sysdat = read_qfile(inpf)
      ress[rk] = bres.copy()
      ress[rk]['kpoint'] = sysdat['system']['kpoint']
    except IOError:
      print "  (cannot find QMC output, skipping)"
      continue
    
    print "  energies..." 
    try:
      inpf = open(kroot+'.dmc.log','r')
      egydat = read_qenergy(inpf,'/home/brian/bin/gosling')
      ress[rk]['dmc_energy']     =  egydat['egy']
      ress[rk]['dmc_energy_err'] =  egydat['err']
    except IOError:
      print "  (cannot find ground state energy)"

    try:
      inpf = open(kroot+'.ogp.log','r')
      ogpdat = read_qenergy(inpf,'/home/brian/bin/gosling')
      ress[rk]['dmc_excited_energy']     =  ogpdat['egy']
      ress[rk]['dmc_excited_energy_err'] =  ogpdat['err']
    except IOError:
      print "  (cannot find excited state energy)"

    print "  fluctuations..." 
    try:
      inpf = open(kroot+'.ppr.o','r')
      fludat, fluerr = read_number_dens(inpf)
      avg, var, cov, avge, vare, cove = moments(fludat,fluerr)
      ress[rk]['covariance']          =  cov
      ress[rk]['site_charge']         =  [] # TODO
    except IOError:
      print "  (cannot find number fluctuation)"

    print "  1-RDM..." 
    try:
      inpf = open(kroot+'.ordm.o')
      odmdat = read_dm(inpf)
      if odmdat==None:
        print "  (Error in 1-RDM output, skipping)"
      else:
        ress[rk]['1rdm'] = odmdat
    except IOError:
      print "  (cannot find 1-RDM)"

    print "  done."; 

  return ress

# Read a directory which has multiple files with data into a single dictionary
# with relevant information.
def read_dir(froot,gosling='./gosling'):
  """ Reads a CRYSTAL + QWalk directory's data into a dictionary.
  
  Current dictionary keys:
    fmixing, kdens, excited_energy_err, ordering, total_energy_err, 
    excited_energy, ts, nfu, tole, total_energy, se_height, tolinteg, 
    supercell, mixing, kpoint, spinlock, dft_energy, a, c, broyden, 
    dft_moments, 1rdm, access_root
  """

  ############################################################################
  # This first section should be edited to reflect naming conventions!
  dftfile   = froot+'.d12'
  dftoutf   = froot+'.d12.out'
  # Start by taking only the real k-points, since I'm sure these are on solid
  # ground, and have enough sample points. TODO generalize
  realk = array([1,3,8,10,27,29,34,36]) - 1
  ############################################################################

  bres = {} # Data that is common to all k-points.
  ress = {} # Dict of all k-point data in directory.

  bres['access_root'] = getcwd() + '/' + froot

  metad = convert_to_metadata(froot)
  bres['ordering'] = metad['mag']


  
  print "Working on",froot+"..." 
  print "  DFT params and results..." 
  try:
    dftdat = read_cryinp(open(dftfile,'r'))
    bres['a'] = dftdat['latparms'][0]
    bres['c'] = dftdat['latparms'][1]
    bres['se_height'] = dftdat['apos'][dftdat['atypes'].index(234)][-1]
    for key in ['mixing','broyden','fmixing','tolinteg',
                'kdens','spinlock','supercell','tole']:
      bres[key] = dftdat[key]
    bres['nfu'] = int(round(linalg.det(array(bres['supercell']).reshape((3,3)))))
    dftdat = read_cryout(open(dftoutf,'r'))
    bres['dft_energy'] = dftdat['dft_energy']
    bres['dft_moments'] = dftdat['dft_moments']
  except IOError:
    print "There's no dft in this directory!"
    return {}

  for rk in realk:
    kroot = froot + '_' + str(rk)

    print "  now DMC:",kroot+"..." 
    try:
      sysdat = read_qfile(open(kroot+'.sys','r'))
      dmcinp = read_qfile(open(kroot+'.dmc','r'))
      ress[rk] = bres.copy()
      ress[rk]['kpoint'] = sysdat['system']['kpoint']
      ress[rk]['ts'] = dmcinp['method']['timestep']
    except IOError:
      print "  (cannot find QMC input, skipping)"
      continue
    
    print "  energies..." 
    try:
      inpf = open(kroot+'.dmc.log','r')
      egydat = read_qenergy(inpf,gosling)
      ress[rk]['dmc_energy']     =  egydat['egy']
      ress[rk]['dmc_energy_err'] =  egydat['err']
    except IOError:
      print "  (cannot find ground state energy log file)"

    try:
      inpf = open(kroot+'.ogp.log','r')
      ogpdat = read_qenergy(inpf,gosling)
      ress[rk]['dmc_excited_energy']     =  ogpdat['egy']
      ress[rk]['dmc_excited_energy_err'] =  ogpdat['err']
    except IOError:
      print "  (cannot find excited state energy log file)"

    print "  fluctuations..." 
    try:
      inpf = open(kroot+'.ppr.o','r')
      fludat, fluerr = read_number_dens(inpf)
      if fludat==None:
        print "  (Error in number fluctuation output, skipping)"
      avg, var, cov, avge, vare, cove = moments(fludat,fluerr)
      ress[rk]['covariance']          =  cov
      ress[rk]['site_charge']         =  [] # TODO
    except IOError:
      print "  (cannot find number fluctuation)"

    print "  1-RDM..." 
    try:
      inpf = open(kroot+'.ordm.o')
      odmdat = read_dm(inpf)
      if odmdat==None:
        print "  (Error in 1-RDM output, skipping)"
      else:
        ress[rk]['1rdm'] = odmdat
    except IOError:
      print "  (cannot find 1-RDM)"

    print "  done."; 

  return ress
