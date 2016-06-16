import numpy as np
import json
import data_processing as dp
import pandas as pd
from pymatgen.io.cif import CifParser

# TODO generalize!
VARTOL = 1e-2
NFE = 8
NORBFE = 10
NORBCH = 4
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
  copykeys = ['dft','supercell','total_spin','charge','xyz','cif','control']
  for copykey in copykeys:
    if copykey in record.keys():
      res[copykey] = record[copykey]
  res['dft'] = record['dft']
  if 'mag_moments' in record['dft'].keys():
    res['dft']['spins_consistent'] = _check_spins(res['dft'],small=SMALLSPIN)
  if 'vmc' in record['qmc'].keys():
    res['vmc'] = _process_vmc(record['qmc']['vmc'])
  res['dmc'] = _process_dmc(record['qmc']['dmc'])
  if 'postprocess' in record['qmc'].keys():
    res['dmc'].update(_process_post(record['qmc']['postprocess']))
    res['dmc'].update(_process_post(record['qmc']['postprocess']))
  return res

def _process_post(post_record):
  """ Process postprocess results by k-averaging and site-averaging."""
  if 'results' not in post_record.keys(): return {}

  res = {}
  # Right now just checks the first k-point: problem?
  if 'region_fluctuation' in post_record['results'][0]['results']['properties'].keys():
    res['fluct'] = _analyze_nfluct(post_record)
  if 'tbdm_basis' in post_record['results'][0]['results']['properties'].keys():
    res['ordm'] = _analyze_ordm(post_record)
  return res

def _process_vmc(dmc_record):
  grouplist = ['jastrow','optimizer']
  res = {}
  if 'results' not in dmc_record.keys():
    return res
  res['energy'] = json.loads(pd.DataFrame(dmc_record['results'])\
      .groupby(grouplist)\
      .apply(_kaverage_energy)\
      .reset_index()
      .to_json()
    )
  return res

def _process_dmc(dmc_record):
  grouplist = ['timestep','jastrow','localization','optimizer']
  res = {}
  if 'results' not in dmc_record.keys():
    return res
  res['energy'] = json.loads(pd.DataFrame(dmc_record['results'])\
      .groupby(grouplist)\
      .apply(_kaverage_energy)\
      .reset_index()
      .to_json()
    )
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
        'variance_err':(sgrp['varerr']**2).mean()**0.5,
        'magmom':abs(sgrp['magmom'].values).mean(),
        'magmom_err':(sgrp['magmom_err']**2).mean()**0.5
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
  avgdf['element'] = "Te"
  avgdf.loc[avgdf['site']<NFE,'element'] = "Fe"

  # Site average.
  savgdf = avgdf.groupby(grouplist+['spinchan','element'])\
      .apply(siteaverage)\
      .reset_index()
  magdf = savgdf.drop(['spinchan','variance','variance_err'],axis=1).drop_duplicates()
  vardf = savgdf.drop(['magmom','magmom_err'],axis=1)
  return { 'magmom':json.loads(magdf.to_json()),
           'variance':json.loads(vardf.to_json()) }

def _analyze_ordm(post_record):
  """ Compute physical values and site-average 1-body RDM. """
  def saverage_orb(sgrp):
    if sgrp['ordm'].std() > VARTOL:
      print("Site average warning: variation in sites larger than expected.")
      print("%.3f > %.2f"%(sgrp['ordm'].std(),VARTOL))
    return pd.Series({
        'ordm':sgrp['ordm'].mean(),
        'ordm_err':(sgrp['ordm_err']**2).mean()**0.5,
      })
  def saverage_hop(sgrp):
    if sgrp['ordm'].std() > VARTOL:
      print("Site average warning: variation in sites larger than expected.")
      print("%.3f > %.2f"%(sgrp['ordm'].std(),VARTOL))
    return pd.Series({
        'ordm':sgrp['ordm'].mean(),
        'ordm_err':(sgrp['ordm_err']**2).mean()**0.5,
      })
  grouplist = ['timestep','jastrow','localization','optimizer']
  # k-average (currently selects gamma-only due to bug).
  ordmdf = pd.DataFrame(post_record['results'])\
      .groupby(['timestep','jastrow','localization','optimizer'])\
      .apply(_kaverage_ordm)\
      .reset_index()
  # Classify orbitals based on index.
  infodf = ordmdf['orbni'].drop_duplicates().apply(lambda orbnum:
      pd.Series(dict(zip(['orbnum','elem','atom','orb'],orbinfo(orbnum)))))
  ordmdf = pd.merge(ordmdf,infodf,how='outer',left_on='orbni',right_on='orbnum')
  ordmdf = pd.merge(ordmdf,infodf,how='outer',left_on='orbnj',right_on='orbnum',
      suffixes=("i","j"))
  ordmdf = ordmdf.drop(['orbnumi','orbnumj'],axis=1)
  # Classify atoms based on spin occupations.
  occdf = ordmdf[ordmdf['orbni']==ordmdf['orbnj']]\
      .groupby(grouplist+['atomi'])\
      .agg({'up':np.sum,'down':np.sum})\
      .reset_index()\
      .rename(columns={'atomi':'at'})
  occdf['net']  = occdf['up'] - occdf['down']
  occdf = occdf.drop(['up','down'],axis=1)
  occdf['atspin'] = 'up'
  occdf.loc[occdf['net'] < 0,'atspin'] = 'down'
  occdf.loc[occdf['net'].abs() < 1e-1,'atspin'] = 'zero'
  ordmdf = pd.merge(ordmdf,occdf,
      left_on=grouplist+['atomi'],right_on=grouplist+['at'])
  ordmdf = pd.merge(ordmdf,occdf,
      left_on=grouplist+['atomj'],right_on=grouplist+['at'],
      suffixes=('i','j'))\
      .drop(['ati','atj'],axis=1)
  ordmdf['rel_atspin'] = "antiparallel"
  ordmdf.loc[ordmdf['atspini']==ordmdf['atspinj'],'rel_atspin'] = "parallel"
  ordmdf.loc[ordmdf['atspini']=='zero','rel_atspin'] = "zero"
  ordmdf.loc[ordmdf['atspinj']=='zero','rel_atspin'] = "zero"
  # Classify spin channels based on minority and majority channels.
  ordmdf = ordmdf.set_index([c for c in ordmdf.columns 
    if c not in ['up','down','up_err','down_err']])
  vals = ordmdf[['up','down']].stack()
  vals.index.names = vals.index.names[:-1]+['spin']
  errs = ordmdf[['up_err','down_err']]\
      .rename(columns={'up_err':'up','down_err':'down'})\
      .stack()
  errs.index.names = errs.index.names[:-1]+['spin']
  ordmdf = pd.DataFrame({'ordm':vals,'ordm_err':errs}).reset_index()
  ordmdf['spini'] = "minority"
  ordmdf['spinj'] = "minority"
  ordmdf.loc[ordmdf['spin'] == ordmdf['atspini'],'spini'] = "majority"
  ordmdf.loc[ordmdf['spin'] == ordmdf['atspinj'],'spinj'] = "majority"
  ordmdf.loc[ordmdf['atspini'] == 'zero','spini'] = 'neither'
  ordmdf.loc[ordmdf['atspinj'] == 'zero','spinj'] = 'neither'
  # Focus in on orbital occupations.
  orboccdf = ordmdf[ordmdf['orbni']==ordmdf['orbnj']]\
      .drop([col for col in ordmdf.columns if col[-1]=='j'],1)\
      .groupby(grouplist+['elemi','orbi','spini'])\
      .apply(saverage_orb)\
      .reset_index()
  # Focus in on parallel or antiparallel hopping.
  orbsumsel = grouplist+['atomi','atomj','elemi','elemj','rel_atspin','spini','spinj']
  siteavgsel = [c for c in orbsumsel if c not in ['atomi','atomj']]
  hopdf = ordmdf[ordmdf['atomi'] != ordmdf['atomj']]\
      .groupby(orbsumsel)\
      .agg({'ordm':lambda x:x.abs().sum(), 'ordm_err':lambda x:sum(x**2)**0.5})\
      .reset_index()\
      .groupby(siteavgsel)\
      .agg({'ordm':np.mean, 'ordm_err':lambda x:np.mean(x**2)**0.5})\
      .reset_index()
  return {'orb':json.loads(orboccdf.to_json()),
          'hop':json.loads(hopdf.to_json())}

def _kaverage_energy(kavgdf):
  # Keep unpacking until reaching energy.
  egydf = \
    unpack(
      unpack(
        unpack(
          kavgdf
        ['results'])\
      ['properties'])\
    ['total_energy']).applymap(dp.unlist)
  # TODO generalize!
  weights = np.tile(1./egydf['value'].shape[0],egydf['value'].shape)
  return pd.Series({
    "value":(weights*egydf['value'].values).sum(),
    "error":((weights*egydf['error'].values)**2).sum()**.5
    })

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
  res = pd.DataFrame(datdf['up'].iloc[0]).stack().to_frame('up')
  res = res.join(pd.DataFrame(datdf['down'].iloc[0]).stack().to_frame('down'))
  res = res.join(pd.DataFrame(datdf['up_err'].iloc[0]).stack().to_frame('up_err'))
  res = res.join(pd.DataFrame(datdf['down_err'].iloc[0]).stack().to_frame('down_err'))
  res = res.reset_index()\
      .rename(columns={'level_0':'orbni','level_1':'orbnj'})\
      .set_index(['orbni','orbnj'])
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
    # Note casting prevents numpy.bool.
    return bool((moms == np.array(init_spins)).all())

def orbinfo(orbnum):
  """ Compute orbital info based on orbital number: [element,atomnum,orbital].

  Currently only useful for Fe-chalcogenides. Warning: this depends on how you
  define the basis!"""
  NFe = 8
  NSe = 8
  # CRYSTAL: 'order of internal storage'.
  # s, px, py, pz, dz2-r2, dxz, dyz, dx2-y2, dxy, ...
  Feorbs = ['3s','3px','3py','3pz','4s','3dz2-r2','3dxz','3dyz','3dx2-y2','3dxy']
  Seorbs = ['3s','3px','3py','3pz']
  NbFe = len(Feorbs)
  NbSe = len(Seorbs)
  res = [orbnum]
  if float(orbnum)/(NFe * NbFe) > (1 - 1e-8):
    res += ['Te',(orbnum - NFe*NbFe) // NbSe + 1 + NFe]
    res.append(Seorbs[orbnum%NbSe])
  else:
    res += ['Fe',orbnum // NbFe + 1]
    res.append(Feorbs[orbnum%NbFe])
  return res

###############################################################################
# Format autogen group of function.
def format_datajson(inp_json="results.json",filterfunc=lambda x:True):
  """ Takes processed autogen json file and organizes it into a Pandas DataFrame."""
  rawdf = pd.read_json(open(inp_json,'r'))
  rawdf['nfu'] = rawdf['supercell'].apply(lambda x:
      2*np.linalg.det(np.array(x))
    )
  # Unpacking the energies.
  dftdf = _format_dftdf(rawdf)
  rawdf = rawdf[dftdf['id'].apply(filterfunc)]
  dmcdf = unpack(rawdf['dmc'])
  if 'energy' in dmcdf.columns:
    dmcdf = dmcdf.join(
          unpack(dmcdf['energy'].dropna()).applymap(dp.undict)
        )
    dmcdf = dmcdf.rename(columns={'value':'dmc_energy','error':'dmc_energy_err'})
  alldf = dmcdf.join(dftdf)
  if 'dmc_energy' in dmcdf.columns:
    alldf['dmc_energy'] = alldf['dmc_energy']/alldf['nfu']
    alldf['dmc_energy_err'] = alldf['dmc_energy_err']/alldf['nfu']
  listcols = [
      'broyden',
      'initial_charges',
      'energy_trace',
      'initial_spin',
      'kmesh'
#      'localization',
#      'timestep',
#      'jastrow',
#      'optimizer'
    ]

  if 'mag_moments' in rawdf.columns: listcols.append('mag_moments')

  # Convert lists.
  for col in listcols:
    alldf.loc[alldf[col].notnull(),col] = \
        alldf.loc[alldf[col].notnull(),col].apply(lambda x:tuple(x))
  for col in alldf.columns:
    alldf[col] = pd.to_numeric(alldf[col],errors='ignore')

  return alldf

def _format_dftdf(rawdf):
  def desect_basis(basis_info):
    if type(basis_info)==list:
      return pd.Series(dict(zip(
        ['basis_lowest','basis_number','basis_factor'],basis_info)))
    elif type(basis_info)==dict:
      min_basis = 1e10
      for atom in basis_info.keys():
        new = min([np.array(elem['coefs'])[0,:].min() for elem in basis_info[atom]])
        if new < min_basis: min_basis = new
      return pd.Series(dict(zip(
        ['basis_lowest','basis_number','basis_factor'],[min_basis,0,0])))
    else:
      return pd.Series(dict(zip(
        ['basis_lowest','basis_number','basis_factor'],[0,0,0])))
  def cast_supercell(sup):
    for rix,row in enumerate(sup):
      sup[rix] = tuple(row)
    return tuple(sup)
  ids = rawdf['control'].apply(lambda x:x['id'])
  dftdf = unpack(rawdf['dft'])
  dftdf = dftdf.join(ids).rename(columns={'control':'id'})
  for rawinfo in ['supercell','nfu','cif','xyz']:
    if rawinfo in rawdf.columns:
      dftdf = dftdf.join(rawdf[rawinfo])
  funcdf = pd.DataFrame(dftdf['functional'].to_dict()).T
  dftdf = dftdf.join(funcdf)
  dftdf['tolinteg'] = dftdf['tolinteg'].apply(lambda x:x[0])
  dftdf['spins_consistent'] = dftdf['spins_consistent'].astype(bool)
  dftdf = dftdf.join(dftdf['basis'].apply(desect_basis))
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
