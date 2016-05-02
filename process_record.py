import numpy as np
import json
import data_processing as dp
import pandas as pd
from pymatgen.io.cif import CifParser

# TODO generalize!
VARTOL = 1e-2
NFE = 8
SMALLSPIN = 1.0 # Spins less than this are considered zero.

################################################################################
# If you're wondering about how to use these, and you're in the Wagner group on
# github, check out my FeTe notebook!
################################################################################

###############################################################################
# Process record group of functions.
def process_record(record):
  """ Take the json produced from autogen and process into a dictionary of much
  processed and more useful results. """
  res = {}
  copykeys = ['dft','supercell','total_spin','charge','cif','control']
  for copykey in copykeys:
    res[copykey] = record[copykey]
  res['dft'] = record['dft']
  if 'mag_moments' in record['dft'].keys():
    res['dft']['spins_consistent'] = _check_spins(res['dft'],small=SMALLSPIN)
  res['dmc'] = _process_dmc(record['qmc']['dmc'])
  res['dmc'].update(_process_post(record['qmc']['postprocess']))
  res['dmc'].update(_process_post(record['qmc']['postprocess']))
  return res

def _process_post(post_record):
  """ Process postprocess results by k-averaging and site-averaging."""
  if 'results' not in post_record.keys(): return {}

  res = {}
  # Right now just checks the first k-point: problem?
  if 'region_fluctuation' in post_record['results'][0]['results']['properties'].keys():
    nfluctdf = _analyze_nfluct(post_record)
    res['fluct'] = json.loads(nfluctdf.reset_index().to_json())
  if 'tbdm_basis' in post_record['results'][0]['results']['properties'].keys():
    ordmdf = _analyze_ordm(post_record)
    res['ordm']  = ordmdf

  return res

def _process_dmc(dmc_record):
  if 'results' not in dmc_record.keys():
    return {}
  res = {}
  for savekey in ['excitations','timestep','jastrow',
      'optimizer','localization','nblock']:
    res[savekey] = dmc_record[savekey]
  res['energy'] = _kaverage_energy(dmc_record['results'])
  return res

def _analyze_nfluct(post_record):
  """ Compute physical values and site-average number fluctuation. """
  def diag_exp(rec):
    """ Compute mean and variance. """
    res = {}
    for dat in ['avg','var','avgerr','varerr']:
      res[dat] = 0.0
    for info in ['jastrow', 'optimizer', 'localization', 
                 'timestep', 'spini', 'sitei']:
      res[info] = rec[info]
    pmat     =  rec['value']
    perr     =  rec['error']
    nmax = len(pmat)
    for n in range(nmax): 
      res['avg']    += n*pmat[n][n]
      res['avgerr'] += (n*perr[n][n])**2
    res['avgerr']= res['avgerr']**0.5
    for n in range(nmax): 
      res['var']    += (n-res['avg'])**2*pmat[n][n]
      res['varerr'] += (perr[n][n]*(n-res['avg'])**2)**2 +\
          (2*pmat[n][n]*res['avgerr']*(n-res['avg']))**2
    res['varerr'] = res['varerr']**0.5
    return pd.Series(res)

  def covar(rec,adf):
    """ Compute covariance. """
    cov = 0.0
    pmat = rec['value']
    nmax = len(pmat)
    avgi = adf.loc[(rec['spini'],rec['sitei']),'avg']
    avgj = adf.loc[(rec['spinj'],rec['sitej']),'avg']
    for m in range(nmax): 
      for n in range(nmax): 
        cov += pmat[m][n]*(m-avgi)*(n-avgj)
    return pd.Series({
        'spini':rec['spini'],
        'spinj':rec['spinj'],
        'sitei':rec['sitei'],
        'sitej':rec['sitej'],
        'cov':cov
      })

  def subspins(siterec):
    tmpdf = siterec.set_index('spin')
    magmom = tmpdf.loc['up','avg'] - tmpdf.loc['down','avg']
    magerr = (tmpdf.loc['up','avgerr']**2 + tmpdf.loc['down','avgerr']**2)**0.5
    return pd.Series({
        'site':siterec['site'].values[0],
        'magmom':magmom, 'magmom_err':magerr
      })

  def siteaverage(sgrp):
    if sgrp['var'].std() > VARTOL:
      print("Site average warning: variation in sites larger than expected.")
      print("%f > 1e-2"%sgrp['var'].std())
    return pd.Series({
        'variance':sgrp['var'].mean(),
        'variance_err':(sgrp['varerr']**1).mean()**0.5,
        'magmom':abs(sgrp['magmom'].values).mean(),
        'magmom_err':(sgrp['magmom']**2).mean()**0.5
      })

  # Moments and other arithmatic.
  #fluctdf = _kaverage_fluct(post_record['results'])
  grouplist = ['timestep','jastrow','localization','optimizer']
  fluctdf = pd.DataFrame(post_record['results'])\
      .groupby(grouplist)\
      .apply(_kaverage_fluct)\
      .reset_index()
  for s in ['spini','spinj']:
    ups = (fluctdf[s] == 0)
    fluctdf[s] = "down"
    fluctdf.loc[ups,s] = "up"
  diag = ( (fluctdf['spini']==fluctdf['spinj']) &\
           (fluctdf['sitei']==fluctdf['sitej'])    )
  avgdf = fluctdf[diag].apply(diag_exp,axis=1)
  avgdf = avgdf.rename(columns={'spini':'spin','sitei':'site'})
  #covdf = fluctdf.apply(lambda x: covar(x,avgdf.set_index(['spin','site'])),axis=1)
  magdf = avgdf.groupby(grouplist+['site']).apply(subspins)
  avgdf = pd.merge(avgdf,magdf)

  # Catagorization.
  avgdf['netmag'] = "down"
  avgdf.loc[avgdf['magmom']>0,'netmag'] = "up"
  avgdf['spinchan'] = "minority"
  avgdf.loc[avgdf['netmag']==avgdf['spin'],'spinchan'] = "majority"
  avgdf['element'] = "Se"
  avgdf.loc[avgdf['site']<NFE,'element'] = "Fe"

  # Site average.
  savgdf = avgdf.groupby(grouplist+['spinchan','element']).apply(siteaverage)

  return savgdf

def _analyze_ordm(post_record):
  """ Compute physical values and site-average 1-body RDM. """
  ordmdf = pd.DataFrame(post_record['results'])\
      .groupby(['timestep','jastrow','localization','optimizer'])\
      .apply(_kaverage_ordm)
  return ordmdf

def _kaverage_energy(reclist):
  # Warning! _kaverage_qmc() assuming equal k-point weight!
  keys = reclist[0].keys()
  # Check if implementation can handle data.
  for option in [k for k in keys if k not in ['results','knum']]:
    for rec in reclist:
      if (type(rec[option])==list) and (len(rec[option]) != 1):
        print("kaverage exception:",rec[option])
        raise AssertionError("Error! _kaverage_qmc() takes no note of timestep or localization lists!"+\
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
    ['total_energy']).applymap(dp.unlist)
  return {"value":egydf['value'].mean(),"error":(egydf['error']**2).mean()**.5}

def _kaverage_fluct(reclist):
  # Warning! _kaverage_qmc() assuming equal k-point weight!
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
  spiniser = datdf.applymap(lambda x: x['spin'][0]).drop_duplicates()
  spinjser = datdf.applymap(lambda x: x['spin'][1]).drop_duplicates()
  siteiser = datdf.applymap(lambda x: x['region'][0]).drop_duplicates()
  sitejser = datdf.applymap(lambda x: x['region'][1]).drop_duplicates()
  valser  = datdf.applymap(lambda x: x['value']).apply(dp.mean_array)
  errser  = datdf.applymap(lambda x: x['error']).apply(dp.mean_array_err)
  # Safely turn DataFrame into Series.
  if spiniser.shape[0] == 1: spiniser = spiniser.iloc[0]
  if spinjser.shape[0] == 1: spinjser = spinjser.iloc[0]
  if siteiser.shape[0] == 1: siteiser = siteiser.iloc[0]
  if sitejser.shape[0] == 1: sitejser = sitejser.iloc[0]
  ret = pd.DataFrame({
      'spini':spiniser,
      'spinj':spinjser,
      'sitei':siteiser,
      'sitej':sitejser,
      'value':valser,
      'error':errser
    }).set_index(['spini','spinj','sitei','sitej'])
  return ret

def _kaverage_ordm(kavgdf):
  # Warning! _kaverage_qmc() assuming equal k-point weight!
  res = {}
  datdf =\
    unpack(
      unpack(
        unpack(
          unpack(
            kavgdf\
          ['results'])\
        ['properties'])\
      ['tbdm_basis'])\
    ['obdm'])
  for spin in ('up','down'):
    #print("datdf[spin]\n",datdf[spin])
    #res[spin] = dp.abs_mean_array(datdf[spin])
    #res[spin+'_err'] = dp.mean_array_err(datdf[spin+'_err'])
    res[spin] = datdf[spin].iloc[4]
    res[spin+'_err'] = datdf[spin+'_err'].iloc[4]
  return res

def _check_spins(dft_record,small=1.0):
  """ Check that the spins that were set at the beginning correspond to the
  spins it ends up having. Anything less than small is considered zero."""
  init_spins = dft_record['initial_spin']
  moms = dft_record['mag_moments']
  moms = np.array(moms)
  zs = abs(moms) < small
  up = moms > 0.
  dn = moms < 0.
  moms.dtype = int
  moms[up] = 1
  moms[dn] = -1
  moms[zs] = 0
  if len(init_spins)==0:
    if (moms == np.zeros(moms.shape)).all():
      return True
    else:
      return False
  else:
    return (moms == np.array(init_spins)).all()

###############################################################################
# Format autogen group of functions (this was used before process_records).
def format_autogen(inp_json="results.json"):
  """
  Takes autogen json file and organizes it into a Pandas DataFrame.
  """
  rawdf = pd.read_json(open(inp_json,'r'))
  rawdf['nfu'] = rawdf['supercell'].apply(lambda x:
      2*np.linalg.det(np.array(x))
    )
  # Unpacking the energies.
  dftdf = _format_dftdf(rawdf)
  dmcdf = unpack(rawdf['dmc'])
  dmcdf = dmcdf.join(unpack(dmcdf['energy']))
  dmcdf = dmcdf.rename(columns={'value':'dmc_energy','error':'dmc_energy_err'})
  alldf = dmcdf.join(dftdf)
  alldf['dmc_energy'] = alldf['dmc_energy']/alldf['nfu']
  alldf['dmc_energy_err'] = alldf['dmc_energy_err']/alldf['nfu']
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

  # Convert lists.
  for col in listcols:
    alldf.loc[alldf[col].notnull(),col] = \
        alldf.loc[alldf[col].notnull(),col].apply(lambda x:tuple(x))
  for col in alldf.columns:
    alldf[col] = pd.to_numeric(alldf[col],errors='ignore')

  # Number fluctuation.
  sel = alldf['fluct'].notnull()
  fluctdf = alldf.loc[sel,'fluct'].apply(lambda df:
    pd.DataFrame(df).set_index(['element','spinchan']).stack())
  alldf = alldf.join(fluctdf)
  alldf = dp.untuple_cols(alldf,"fluct")

  return alldf

def _format_dftdf(rawdf):
  def desect_basis(df):
    return pd.Series(dict(zip(['basis_lowest','basis_number','basis_factor'],df)))
  def cast_supercell(sup):
    for rix,row in enumerate(sup):
      sup[rix] = tuple(row)
    return tuple(sup)
  ids = rawdf['control'].apply(lambda x:x['id'])
  dftdf = unpack(rawdf['dft'])
  dftdf = dftdf.join(ids).rename(columns={'control':'id'})
  for rawinfo in ['supercell','nfu','cif']:
    dftdf = dftdf.join(rawdf[rawinfo])
  funcdf = pd.DataFrame(dftdf['functional'].to_dict()).T
  dftdf = dftdf.join(funcdf)
  dftdf['tolinteg'] = dftdf['tolinteg'].apply(lambda x:x[0])
  dftdf['spins_consistent'] = dftdf['spins_consistent'].astype(bool)
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

###############################################################################
# Misc. tools.
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
