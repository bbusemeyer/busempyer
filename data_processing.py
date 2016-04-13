#!/usr/bin/python
import numpy as np
import pandas as pd
import os
import json
import cryfiles_io as cio
import qfiles_io as qio
import qefiles_io as qeio
import cubetools as ct
from pymatgen.io.cif import CifParser

################################################################################
# If you're wondering about how to use these, and you're in the Wagner group on
# github, check out my FeTe notebook!
################################################################################

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
    return tuple((abs(ser.apply(np.array)).sum()/ser.shape[0]).tolist())

def mean_array_err(ser):
  """
  Computes the error corresponding to mean_array().
  """
  if ser.isnull().any():
    return None
  else:
    return tuple(((ser.apply(lambda x: (np.array(x)/ser.shape[0])**2)).sum()**.5).tolist())

# After processing, what are the columns that define the accuracy level of the
# calculation? Defined as a function to make it more readable when importing.
# Pls keep alphabetize (<range> sort u). Pls. 
def gen_autogen_defining_columns():
  return [
      'basis_factor',
      'basis_lowest',
      'basis_number',
      'a','b','c',
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

# Temporary function to convert file names to metadata about the calculations.
# In the future, an optional metadata file should contain overall qualitiative
# descriptions of the data, like magnetic ordering, or pressure. These things
# are redundant, but convenient for interpreting sets of numbers together.
def convert_to_metadata(froot):
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

# Read a directory which has multiple files with data into a single dictionary
# with relevant information.
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
    dftdat = cio.read_cryinp(open(dftfile,'r'))
    bres['a'] = dftdat['latparms'][0]
    bres['c'] = dftdat['latparms'][1]
    bres['se_height'] = dftdat['apos'][dftdat['atypes'].index(234)][-1]
    for key in ['mixing','broyden','fmixing','tolinteg',
                'kdens','spinlock','supercell','tole','basis']:
      bres[key] = dftdat[key]
    dftdat = cio.read_cryout(open(dftoutf,'r'))
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
      sysdat = qio.read_qfile(open(kroot+'.sys','r'))
      dmcinp = qio.read_qfile(open(kroot+'.dmc','r'))
      ress[rk] = bres.copy()
      ress[rk]['kpoint'] = sysdat['system']['kpoint']
      ress[rk]['ts'] = dmcinp['method']['timestep']
    except IOError:
      print("  (cannot find QMC input, skipping)")
      continue
    
    print("  energies..." )
    try:
      inpf = open(kroot+'.dmc.log','r')
      egydat = qio.read_qenergy(inpf,gosling)
      ress[rk]['dmc_energy']     =  egydat['egy']
      ress[rk]['dmc_energy_err'] =  egydat['err']
    except IOError:
      print("  (cannot find ground state energy log file)")

    try:
      inpf = open(kroot+'.ogp.log','r')
      ogpdat = qio.read_qenergy(inpf,gosling)
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
      fludat, fluerr = qio.read_number_dens(inpf)
      if fludat is None:
        print("  (Error in number fluctuation output, skipping)")
      else:
        avg, var, cov, avge, vare, cove = qio.moments(fludat,fluerr)
        ress[rk]['average']    = avg
        ress[rk]['covariance'] = cov
    except IOError:
      print("  (cannot find number fluctuation)")

    print("  1-RDM..." )
    try:
      inpf = open(kroot+'.ordm.o')
      odmdat = qio.read_dm(inpf)
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

# Helper for format_autogen().
def _format_dftdf(rawdf):
  def desect_basis(df):
    return pd.Series(dict(zip(['basis_lowest','basis_number','basis_factor'],df)))
  def cast_supercell(sup):
    for rix,row in enumerate(sup):
      sup[rix] = tuple(row)
    return tuple(sup)
  ids = rawdf['control'].apply(lambda x:x['id'])
  dftdf = pd.DataFrame(rawdf['dft'].to_dict()).T
  dftdf = dftdf.join(ids).rename(columns={'control':'id'})
  for rawinfo in ['supercell','nfu','cif']:
    dftdf = dftdf.join(rawdf[rawinfo])
  funcdf = pd.DataFrame(dftdf['functional'].to_dict()).T
  dftdf = dftdf.join(funcdf)
  dftdf['tolinteg'] = dftdf['tolinteg'].apply(lambda x:x[0])
  dftdf = dftdf.join(dftdf['basis'].apply(desect_basis))
  dftdf['basis_number'] = dftdf['basis_number'].apply(lambda x:int(round(x)))
  dftdf['basis_factor'] = dftdf['basis_factor'].apply(lambda x:int(round(x)))
  dftdf.loc[dftdf['supercell'].notnull(),'supercell'] = \
      dftdf.loc[dftdf['supercell'].notnull(),'supercell']\
      .apply(lambda x:cast_supercell(x))
  if 'mag_moments' in dftdf.columns:
    dftdf['max_mag_moment'] = np.nan
    dftdf.loc[dftdf['mag_moments'].notnull(),'max_mag_moment'] =\
        dftdf.loc[dftdf['mag_moments'].notnull(),'mag_moments'].apply(lambda x:
            max(abs(np.array(x)))
          )
  dftdf['dft_energy'] = dftdf['total_energy']/dftdf['nfu']
  for redundant in ['basis','functional']:
    del dftdf[redundant]
  return dftdf

def unpack(ser):
  """ Attempt to turn a series of dictionaries into a new DataFrame. 
  Works with most autogen levels of data storage. """
  return pd.DataFrame(ser.to_dict()).T

# Tuple to DF entry.
def parse_err(df,key='energy'):
  tmpdf = df[key].apply(lambda x: pd.Series({key:x[0],'energy_err':x[1]}))
  del df[key]
  return df.join(tmpdf)

# Get out atomic positions (and possibly more later).
# Pretty slow: can be made faster by saving cifs that are already done.
def extract_struct(cifstr):
  parser = CifParser.from_string(cifstr)\
      .get_structures()[0]\
      .as_dict()
  lat_a = parser['lattice']['a']
  lat_b = parser['lattice']['b']
  lat_c = parser['lattice']['c']
  poss = [
      tuple(site['abc']) for site in 
      parser['sites']
    ]
  ions = [
      site['species'][0]['element'] for site in 
      parser['sites']
    ]
  positions = {}
  for iidx,ion in enumerate(ions):
    if ion in positions.keys():
      positions[ion].append(poss[iidx])
    else:
      positions[ion] = [poss[iidx]]
  for key in positions.keys():
    positions[key] = np.array(positions[key])
  return pd.Series(
      [lat_a,lat_b,lat_c,positions],
      ['a','b','c','positions']
    )

def format_autogen(inp_json="results.json"):
  """
  Takes autogen json file and organizes it into a Pandas DataFrame.
  """
  rawdf = pd.read_json(open(inp_json,'r'))
  rawdf['nfu'] = rawdf['supercell'].apply(lambda x:
      2*np.linalg.det(np.array(x))
    )
  dftdf = _format_dftdf(rawdf)
  qmcdf = unpack(rawdf['qmc'])
  dmcdf = unpack(qmcdf['dmc'])
  alldf = dmcdf.join(dftdf)
  if 'vmc' in qmcdf.columns:
    vmcdf = unpack(qmcdf['vmc'])
    alldf = alldf.join(vmcdf,rsuffix="_vmc")
  if 'postprocess' in qmcdf.columns:
    postdf = unpack(qmcdf['postprocess'])
    alldf = alldf.join(postdf,rsuffix="_post")
  listcols = [
      'broyden',
      'initial_charges',
      'initial_spin',
      'kmesh',
      'localization',
      'timestep',
      'jastrow',
      'optimizer'
    ]
  if 'mag_moments' in rawdf.columns: listcols.append('mag_moments')
  for col in listcols:
    alldf.loc[alldf[col].notnull(),col] = \
        alldf.loc[alldf[col].notnull(),col].apply(lambda x:tuple(x))

  for col in alldf.columns:
    alldf[col] = pd.to_numeric(alldf[col],errors='ignore')

  return alldf

# Safely take list of one element into it's value.
def unlist(li):
  if len(li) > 1: AssertionError("unlist can't operate on multi-element list")
  return li[0]

def kaverage_energy(reclist):
  print("Warning! kaverage_qmc() assuming equal k-point weight!")
  keys = reclist[0].keys()
  # Check if implementation can handle data.
  for option in [k for k in keys if k not in ['results','knum']]:
    for rec in reclist:
      if (type(rec[option])==list) and (len(rec[option]) != 1):
        print(rec[option])
        AssertionError("Error! kaverage_qmc() takes no note of timestep or localization lists!"+\
          "If you need this functionality, I encourage you to generize it for me"+\
          "(and anyone else using it)!")
  # Keep unpacking until reaching energy.
  egydf = \
    unpack(
      unpack(
        unpack(
          pd.DataFrame(reclist)\
        ['results'])\
      ['properties'])\
    ['total_energy']).applymap(unlist)
  return {"value":egydf['value'].mean(),"error":(egydf['error']**2).mean()**.5}

def kaverage_fluct(reclist):
  print("Warning! kaverage_qmc() assuming equal k-point weight!")
  keys = reclist[0].keys()
  # Check if implementation can handle data.
  for option in [k for k in keys if k not in ['results','knum']]:
    for rec in reclist:
      if (type(rec[option])==list) and (len(rec[option]) != 1):
        print(rec[option])
        AssertionError("Error! kaverage_qmc() takes no note of timestep or localization lists!"+\
          "If you need this functionality, I encourage you to generize it for me"+\
          "(and anyone else using it)!")
  # Keep unpacking until reaching energy.
  datdf = \
    unpack(
      unpack(
        unpack(
          unpack(
            pd.DataFrame(reclist)\
          ['results'])\
        ['properties'])\
      ['region_fluctuation'])\
    ['fluctuation data'])
  spinser = datdf.applymap(lambda x: tuple(x['spin'])).drop_duplicates()
  siteser = datdf.applymap(lambda x: tuple(x['region'])).drop_duplicates()
  valser  = datdf.applymap(lambda x: x['value']).apply(mean_array)
  errser  = datdf.applymap(lambda x: x['error']).apply(mean_array_err)
  if spinser.shape[0] == 1: spinser = spinser.T[0]
  if siteser.shape[0] == 1: siteser = siteser.T[0]
  fluctdf = pd.DataFrame([spinser,siteser,valser,errser],
                        ["spin", "site", "value", "error"])
  return fluctdf.T.to_dict()

def process_dmc(dmc_record):
  res = {}
  res['energy'] = kaverage_energy(dmc_record['results'])
  return res

def process_post(post_record):
  res = {}
  res['fluct'] = kaverage_fluct(post_record['results'])
  return res


# Convert pandas DataFrame row into a dictionary.
def row_to_dict(row):
  ret = row.T.to_dict()
  return ret[list(ret.keys())[0]]

