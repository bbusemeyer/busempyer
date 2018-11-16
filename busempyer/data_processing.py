import numpy as np
import os
import json
from busempyer.cryfiles_io import read_cryinp,read_cryout
from busempyer.qfiles_io import read_qfile,read_qenergy,read_number_dens,moments,read_dm
#import read_numberfluct as rn
from busempyer.cubetools import read_cube
from copy import deepcopy

def mean_err(x):
  return sum(x**2)**0.5

# A lot of this is probably obsolete with latest version of autogen.

###############################################################################
# Misc. tools.
def mean_array(ser):
  """
  Converts an interable to an array, averages the series, and returns the array
  of averaged elements.
  """
  if ser.isnull().any():
    return None
  else:
    return (ser.apply(np.array).sum()/ser.shape[0]).tolist()

def abs_mean_array(ser):
  """
  Converts an interable to an array, absolute value averages the series, and
  returns the array of averaged elements.
  """
  if ser.isnull().any():
    return None
  else:
    #return tuple(map(tuple,ser.apply(np.array).sum()))
    return (abs(ser.apply(np.array)).sum()/ser.shape[0]).tolist()

def mean_array_err(ser):
  """
  Computes the error corresponding to mean_array().
  """
  if ser.isnull().any():
    return None
  else:
    return ((ser.apply(lambda x: (np.array(x)/ser.shape[0])**2)).sum()**.5).tolist()

def gen_autogen_defining_columns():
  """ After processing, what are the columns that define the accuracy level of the
  calculation? Defined as a function to make it more readable when importing.
  Pls keep alphabetize (<range> sort u).  """
  return [
      'basis_factor',
      'basis_number',
      'basis_lowest',
      'levshift',
      'a','b','c',
      'total_spin',
      'correlation',
      'exchange',
      'hybrid',
      'jastrow',
      'kmesh',
      'localization',
      'magstate',
      'optimizer',
      'struct',
      'supercell',
      'timestep',
      'tolinteg'
    ]

def convert_to_metadata(froot):
  """ Temporary function to convert file names to metadata about the calculations.
  In the future, an optional metadata file should contain overall qualitiative
  descriptions of the data, like magnetic ordering, or pressure. These things
  are redundant, but convenient for interpreting sets of numbers together.  """
  mconv = {'che':'checkerboard',
           'str':'collinear',
           'unp':'unpolarized',
           'fer':'ferromagnetic',
           'bic':'bicollinear',
           'fst':'collinear, flip 1',
           'dim':'collinar, flip 2',
           'sta':'staggered',
           'cpx':'checkerboard',
           'spx':'collinear',
           'bpx':'bicollinear'}
  basename = froot.split('/')[-1]
  try:
    mag = mconv[basename[:3]]
  except:
    mag = None
  try:
    if 'px' in basename[:3]: # e.g. 'cpx'
      prs = float(basename.split('_')[1])
    else:
      prs = float(basename[3:].replace('_prs',''))
  except:
    prs = None
  return {'magnetic ordering':mag,'pressure':prs}

def read_dir(froot,gosling='./gosling',read_cubes=False):
  """ Reads a CRYSTAL + QWalk directory's data into a dictionary.
  
  Current dictionary keys:
    fmixing, kdens, excited_energy_err, ordering, total_energy_err, 
    excited_energy, ts, tole, total_energy, se_height, tolinteg, 
    supercell, mixing, kpoint, spinlock, dft_energy, a, c, broyden, 
    dft_moments, 1rdm, access_root, average, covariance
  """

  ############################################################################
  # This first section should be edited to reflect naming conventions!
  dftfile   = froot+'.d12'
  dftoutf   = froot+'.d12.out'
  # Start by taking only the real k-points, since I'm sure these are on solid
  # ground, and have enough sample points. TODO generalize
  realk1 = np.array([1,3,8,10,27,29,34,36]) - 1
  realk2 = np.array([1,4,17,20,93,96,109,112]) - 1
  oldrealk = np.array(['k0','k1','k2','k3','k4','k5','k6','k7'])
  ############################################################################

  bres = {} # Data that is common to all k-points.
  ress = {} # Dict of all k-point data in directory.

  bres['access_root'] = os.getcwd() + '/' + froot

  print("Working on",froot+"..." )
  try:
    with open(froot+'_metadata.json','r') as metaf:
      metad = json.load(metaf)
    try:
      if metad['magnetic ordering'] != None:
        bres['ordering'] = metad['magnetic ordering']
    except KeyError: pass
    try:
      if metad['pressure'] != None:
        bres['pressure'] = metad['pressure']
    except KeyError: pass
    try:
      if metad['pressure'] != None:
        bres['pressure'] = metad['pressure']
    except KeyError: pass
  except IOError:
    print("  Didn't find any metadata")
  
  print("  DFT params and results..." )
  try:
    dftdat = read_cryinp(open(dftfile,'r'))
    bres['a'] = dftdat['latparms'][0]
    bres['c'] = dftdat['latparms'][1]
    bres['se_height'] = dftdat['apos'][dftdat['atypes'].index(234)][-1]
    for key in ['mixing','broyden','fmixing','tolinteg',
                'kdens','spinlock','supercell','tole','basis']:
      bres[key] = dftdat[key]
    dftdat = read_cryout(open(dftoutf,'r'))
    bres['dft_energy'] = dftdat['dft_energy']
    bres['dft_moments'] = dftdat['dft_moments']
  except IOError:
    print("There's no dft in this directory!")
    return {}

  # Determine k-point set and naming convention.
  if os.path.isfile(froot+'_'+str(oldrealk[0])+'.sys'):
    realk=oldrealk
    print("Using old-style kpoint notation.")
  elif os.path.isfile(froot+'_'+str(realk2[-1])+'.sys'):
    realk=realk2
    print("Using 6x6x6 kpoint notation.")
  else:
    realk=realk1
    print("Using 4x4x4 kpoint notation.")

  for rk in realk:
    kroot = froot + '_' + str(rk)

    print("  now DMC:",kroot+"..." )
    try:
      sysdat = read_qfile(open(kroot+'.sys','r'))
      dmcinp = read_qfile(open(kroot+'.dmc','r'))
      ress[rk] = bres.copy()
      ress[rk]['kpoint'] = sysdat['system']['kpoint']
      ress[rk]['ts'] = dmcinp['method']['timestep']
    except IOError:
      print("  (cannot find QMC input, skipping)")
      continue
    
    print("  energies..." )
    try:
      inpf = open(kroot+'.dmc.log','r')
      egydat = read_qenergy(inpf,gosling)
      ress[rk]['dmc_energy']     =  egydat['egy']
      ress[rk]['dmc_energy_err'] =  egydat['err']
    except IOError:
      print("  (cannot find ground state energy log file)")

    try:
      inpf = open(kroot+'.ogp.log','r')
      ogpdat = read_qenergy(inpf,gosling)
      ress[rk]['dmc_excited_energy']     =  ogpdat['egy']
      ress[rk]['dmc_excited_energy_err'] =  ogpdat['err']
    except IOError:
      print("  (cannot find excited state energy log file)")

    if read_cubes:
      print("  densities..." )
      try:
        inpf = open(kroot+'.dmc.up.cube','r')
        upcube = ct.read_cube(inpf)
        inpf = open(kroot+'.dmc.dn.cube','r')
        dncube = ct.read_cube(inpf)
        ress[rk]['updens'] = upcube
        ress[rk]['dndens'] = dncube
      except IOError:
        print("  (cannot find electron density)")
      except ValueError:
        print("  (electron density is corrupted)")

    print("  fluctuations..." )
    try:
      inpf = open(kroot+'.ppr.o','r')
      fludat, fluerr = read_number_dens(inpf)
      if fludat is None:
        print("  (Error in number fluctuation output, skipping)")
      else:
        avg, var, cov, avge, vare, cove = moments(fludat,fluerr)
        ress[rk]['average']    = avg
        ress[rk]['covariance'] = cov
    except IOError:
      print("  (cannot find number fluctuation)")

    print("  1-RDM..." )
    try:
      inpf = open(kroot+'.ordm.o')
      odmdat = read_dm(inpf)
      if odmdat is None:
        print("  (Error in 1-RDM output, skipping)")
      else:
        ress[rk]['1rdm'] = odmdat
    except IOError:
      print("  (cannot find 1-RDM)")

    print("  done." )

  if ress == {}:
    ress['dft-only'] = bres
  return ress

def read_dir_autogen(froot,gosling='./gosling',read_cubes=False):
  """ Reads a CRYSTAL + QWalk directory's data into a dictionary, formats it
  like how autogen would do. NOTE: this is *not* intended to read an
  autogen-generated directory! This is for the purpose of joining an
  autogen-generated database with another set of DMC data.
  
  Current dictionary keys:
    fmixing, kdens, excited_energy_err, ordering, total_energy_err, excited_energy, ts, tole, total_energy, se_height, tolinteg, 
    supercell, mixing, kpoint, spinlock, dft_energy, a, c, broyden, 
    dft_moments, 1rdm, access_root, average, covariance
  """

  ############################################################################
  # This first section should be edited to reflect naming conventions!
  dftfile   = froot+'.d12'
  dftoutf   = froot+'.d12.out'
  # Start by taking only the real k-points, since I'm sure these are on solid
  # ground, and have enough sample points. TODO generalize
  realk1 = list(map(str,np.array([1,3,8,10,27,29,34,36]) - 1))
  realk2 = list(map(str,np.array([1,4,17,20,93,96,109,112]) - 1))
  oldrealk = ['k0','k1','k2','k3','k4','k5','k6','k7']
  ############################################################################

  res  = {}
  res['dft'] = {}
  res['qmc'] = {}

  res['access_root'] = os.getcwd() + '/' + froot

  print("Working on",froot+"..." )
  try:
    with open(froot+'_metadata.json','r') as metaf:
      metad = json.load(metaf)
    try:
      if metad['magnetic ordering'] != None:
        res['ordering'] = metad['magnetic ordering']
    except KeyError: pass
    try:
      if metad['pressure'] != None:
        res['pressure'] = metad['pressure']
    except KeyError: pass
  except IOError:
    print("  Didn't find any metadata")
  
  print("  DFT params and results..." )
  try:
    dftdat = cio.read_cryinp(open(dftfile,'r'))
    res['a'] = dftdat['latparms'][-1]
    res['c'] = dftdat['latparms'][1]
    res['se_height'] = dftdat['apos'][dftdat['atypes'].index(234)][-1]
    res['supercell'] = dftdat['supercell']
    res['total_spin'] = dftdat['spinlock']
    #res['initial_spin'] = dftdat['atomspin']
    res['charge'] = 0
    res['cif'] = "None"
    res['control'] = {'id':froot}
    for key in ['mixing','broyden','fmixing','tolinteg',
                'kdens','tole','basis','initial_spin']:
      res['dft'][key] = dftdat[key]
    dftdat = cio.read_cryout(open(dftoutf,'r'))
    res['dft']['energy'] = dftdat['dft_energy']
    res['dft']['mag_moments'] = dftdat['dft_moments']
  except IOError:
    print("There's no dft in this directory!")
    return res

  # Determine k-point set and naming convention.
  if os.path.isfile(froot+'_'+str(oldrealk[0])+'.sys'):
    realk=oldrealk
    print("Using old-style kpoint notation.")
  elif os.path.isfile(froot+'_'+str(realk2[-1])+'.sys'):
    realk=realk2
    print("Using 6x6x6 kpoint notation.")
  else:
    realk=realk1
    print("Using 4x4x4 kpoint notation.")

  dmc_ret = []
  ppr_ret = []
  for rk in realk:
    entry = {}
    entry['knum'] = int(rk.replace("k",""))
    kroot = froot + '_' + str(rk)

    print("  now DMC:",kroot+"..." )
    try:
      sysdat = read_qfile(open(kroot+'.sys','r'))
      dmcinp = read_qfile(open(kroot+'.dmc','r'))
      entry['kpoint'] = sysdat['system']['kpoint']
      entry['timestep'] = dmcinp['method']['timestep']
      entry['localization'] = "None"
      entry['jastrow'] = "twobody"
      entry['optimizer'] = "variance"
      entry['excitations'] = "no"
    except IOError:
      print("  (cannot find QMC input, skipping)")
      continue
    
    print("  DMC results..." )
    try:
      os.system("{0} -json {1}.dmc.log > {1}.json".format(gosling,kroot))
      dmc_entry = deepcopy(entry)
      dmc_entry['results'] = json.load(open("{}.json".format(kroot)))
      dmc_ret.append(dmc_entry)
    except ValueError:
      print("  (cannot find ground state energy log file)")

    # I'm going to skip the optical gap information for now. Is this really
    # taken care of in autogen yet?
    #try:
    #  inpf = open(kroot+'.ogp.log','r')
    #  ogpdat = qio.read_qenergy(inpf,gosling)
    #  ress[rk]['dmc_excited_energy']     =  ogpdat['egy']
    #  ress[rk]['dmc_excited_energy_err'] =  ogpdat['err']
    #except IOError:
    #  print("  (cannot find excited state energy log file)")

    # I'm going to skip this for now since this is too large for the database
    # generally.
    #if read_cubes:
    #  print("  densities..." )
    #  try:
    #    inpf = open(kroot+'.dmc.up.cube','r')
    #    upcube = ct.read_cube(inpf)
    #    inpf = open(kroot+'.dmc.dn.cube','r')
    #    dncube = ct.read_cube(inpf)
    #    ress[rk]['updens'] = upcube
    #    ress[rk]['dndens'] = dncube
    #  except IOError:
    #    print("  (cannot find electron density)")
    #  except ValueError:
    #    print("  (electron density is corrupted)")

    print("  Postprocessing results..." )
    try:
      inpf = open(kroot+'.ppr.o','r')
      # TODO
      # Apparently this directory is now missing, so you should fix this if
      # needed.
      #fludat = rn.read_number_dens_likejson(inpf)
      if fludat is None:
        print("  (Error in number fluctuation output, skipping)")
      else:
        ppr_entry = deepcopy(entry)
        ppr_entry['results'] = {
            'properties':{
              'region_fluctuation':{
                'fluctuation data': fludat} 
              }
            }
        ppr_ret.append(ppr_entry)
    except IOError:
      print("  (cannot find number fluctuation)")

    # Going to skip for now, since it's not working for my format_autogen yet.
    #print("  1-RDM..." )
    #try:
    #  inpf = open(kroot+'.ordm.o')
    #  odmdat = qio.read_dm(inpf)
    #  if odmdat is None:
    #    print("  (Error in 1-RDM output, skipping)")
    #  else:
    #    ress[rk]['1rdm'] = odmdat
    #except IOError:
    #  print("  (cannot find 1-RDM)")

    print("  done." )

  if len(dmc_ret) > 0:
    res['qmc']['dmc'] = {}
    res['qmc']['dmc']['results'] = dmc_ret
    res['qmc']['dmc']['kpoint'] = sysdat['system']['kpoint']
    res['qmc']['dmc']['timestep'] = dmcinp['method']['timestep']
    res['qmc']['dmc']['localization'] = "None"
    res['qmc']['dmc']['jastrow'] = "twobody"
    res['qmc']['dmc']['optimizer'] = "variance"
    res['qmc']['dmc']['nblock'] = -1
    res['qmc']['dmc']['excitations'] = "no"
  if len(ppr_ret) > 0:
    res['qmc']['postprocess'] = {}
    res['qmc']['postprocess']['results'] = ppr_ret
    res['qmc']['postprocess']['kpoint'] = sysdat['system']['kpoint']
    res['qmc']['postprocess']['timestep'] = dmcinp['method']['timestep']
    res['qmc']['postprocess']['localization'] = "None"
    res['qmc']['postprocess']['jastrow'] = "twobody"
    res['qmc']['postprocess']['optimizer'] = "variance"
    res['qmc']['postprocess']['excitations'] = "no"
    res['qmc']['postprocess']['nblock'] = -1
  return res

def trace_analysis(dftfns,ids=[]):
  """
  Generate a dictionary of energy traces from DFT output, useful for checking
  how certain parameters affect convergence.
  """
  res = {}
  if ids == []: ids = dftfns
  for di,dftfn in enumerate(dftfns):
    with open(dftfn,'r') as dftf:
      res[ids[di]] = cio.read_crytrace(dftf)
  return res

def untuple_cols(df,prefix="",sep="_"):
  replace = {}
  for col in list(df.columns):
    if type(col) == tuple:
      replace[col] = prefix+sep+sep.join(col)
  return df.rename(columns=replace)

# Safely take list of one element into it's value.
def unlist(li):
  if li != li: return np.nan
  assert not len(li) > 1,"unlist can't operate on multi-element list"
  return li[0]

# Safely take list of one element into it's value.
def undict(di):
  if di != di: return np.nan
  assert not len(di) > 1,"undict can't operate on multi-element list"
  return list(di.items())[0][1]

# Convert pandas DataFrame row into a dictionary.
def row_to_dict(row):
  ret = row.T.to_dict()
  return ret[list(ret.keys())[0]]

def cif_to_dict(cifstr):
  """
  Takes a cif string and parses it into a dictionary.

  Currently acts as a inefficient, but effective wrapper around pymatgen's cif
  file parser.
  """
  with open("tmp",'w') as outf:
    outf.write(cifstr)
  cifp = CifParser("tmp").as_dict()
  os.remove("tmp")
  return cifp[cifp.keys()[0]]

def load_duplicate_method_keys(key,jsonstr):
  """ Go through a gosling JSON output, and correct if there are multiple
  instances of this key. Very hacky."""
  json_chunks = jsonstr.split('"')
  key_count = 0
  for chx,chunk in enumerate(json_chunks):
    if chunk == key:
      json_chunks[chx] = "%s_%d"%(key,key_count)
      key_count+=1
  jsonstr = '"'.join(json_chunks).replace('\n','')
  rec = json.loads(jsonstr)
  return rec
###############################################################################

