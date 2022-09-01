''' Misc tools for interfacing with AFQMCLab.'''
import os
from pandas import DataFrame
from busempyer.qmcdata import estimate_warmup
from h5py import File
from numpy import arange
from qharv.reel.forlib.stats import corr as compute_autocorr

def main():
  # Test Paul's autocorrelation calculation. 
  data, safe = read_raw_afqmc("testdata/dt0.010_fc1.00_lec2.00_s0")
  import sys
  sys.path.append("/mnt/home/bbusemeyer/soft/harvest_qmcpack")
  acortime = compute_autocorr(data['energy'].values)

  from pyblock.blocking import find_optimal_block, reblock
  # Compare to reblocking.
  blocked = reblock(data['energy'].values)
  print(blocked)
  optimal = find_optimal_block(len(data),blocked)[0]
  blockedtime = len(data)//blocked[optimal].ndata

  print(f"Autocorrlation time: {acortime}\nBlocked time: {blockedtime}")

def read_afqmc_param(fn):
  afparams = {}
  for line in open(fn,'r').read().split('\n'):
    words = line.split()
    if len(words) != 2: continue
    afparams[words[0]] = words[1]

  return afparams

def dump_afqmc_param(fn='./afqmc_param', **opts):
  ''' dump opts to afqmc_param for AFQMCLab.'''
  outlines = [f"{key} {opts[key]}" for key in opts]
  with open(fn,'w') as outf:
    outf.write('\n'.join(outlines))

def read_afqmc(loc='./',warmup=None,return_trace=False):
  ''' Read AFQMC results in this directory. 

  Args:
    warmup (int): Throw out this many observations before analysis. Default is to estimate based on final half of observations. 
    return_trace (bool): Whether you want the full trace in the results (can be large). 
  Returns:
    results (dict): Block-averaged means and estimates of error.
  '''
  if loc[-1]!='/': loc+='/'

  if not os.path.exists(f"{loc}HNum.dat"):
    print("No results.")
    return {}

  if not os.path.exists(f"{loc}beta.dat"):
    return read_measure_afqmc(loc)

  edf,safe_energy = read_raw_afqmc(loc)

  if warmup is None:
    warmup = estimate_warmup(edf['energy'].values)
    if warmup < 0:
      print(f"\nInsufficient data detected in {loc+'HNum.dat'}.")
      return {
          'warmup': -1,
          'safe_energy': safe_energy,
          'energy': None,
          'stdev': None,
          'error': None,
          'autocortime': 0.0,
          #'blockbeta': None,
          #'blockdata': None
        }


  warmdf = edf.iloc[warmup:]
  autotime = compute_autocorr(warmdf['energy'].values)

  results =  {
      'warmup':       warmup,
      'safe_energy':  safe_energy,
      'energy':       warmdf['energy'].mean(),
      'stdev':        warmdf['energy'].std(),
      'autocortime':  autotime,
      #'blockbeta':    (edf.iloc[warmup:]['beta'].values[blocknbeta//2::blocknbeta]).tolist(),
      #'blockdata':    blockdata['value'].values.tolist(),
    }
  results['error'] = results['stdev']*(autotime/warmdf.shape[0])**0.5

  if return_trace:
    results['beta']  = edf['beta'].values.tolist()
    results['trace'] = edf['energy'].values.tolist()

  return results 

def test_blocking(loc="./",warmup=None):
  ''' Check blocking routine in various ways.'''
  from json import dumps
  results = read_afqmc(loc,warmup)
  print(dumps(results,indent='  '))

  edf,safe_energy = read_raw_afqmc(loc)
  edf = edf.iloc[results['warmup']:]

  from pyqmc.reblock import optimally_reblocked
  pyqmc_results = optimally_reblocked(edf[['energy']])
  print(pyqmc_results)

def read_raw_afqmc(loc="./"):
  if loc[-1] != '/': loc += '/'
  edf = DataFrame([l.split() for l in open(f"{loc}HNum.dat",'r').readlines()],
      columns=('energy','imenergy'),dtype=float)
  edf = edf.join(DataFrame([l.split() for l in open(f"{loc}den.dat",'r').readlines()],
      columns=('weight','imweight'),dtype=float))
  edf = edf.join(DataFrame([l.split() for l in open(f"{loc}beta.dat",'r').readlines()],
      columns=['beta'],dtype=float))
  edf['beta'] = edf['beta'].round(10)
  if edf['beta'].shape != edf['beta'].unique().shape: 
    print("\nWARN (read_raw_afqmc): Repeating betas, hopefully this is a restart. \nFixing betas to be monotonic.")
    dbeta = edf.at[1,'beta'] - edf.at[0,'beta']
    edf['beta'] = arange(edf.shape[0])*dbeta + edf.at[0,'beta']

  # Note: these lines requires population to be large enough:
  safe_energy = edf['energy'].sum() / edf['weight'].sum()
  edf['energy'] = edf['energy'] / edf['weight']

  return edf,safe_energy

def read_measure_afqmc(loc='./'):
  ''' Same as read_afqmc but for measurement without projection. '''
  edf = DataFrame([l.split() for l in open(f"{loc}HNum.dat",'r').readlines()],
      columns=('energy','imenergy'),dtype=float)
  edf = edf.join(DataFrame([l.split() for l in open(f"{loc}den.dat",'r').readlines()],
      columns=('weight','imweight'),dtype=float))
  energy = float([l.split() for l in open(f"{loc}HNum.dat",'r').readlines()][0][0])
  weight = float([l.split() for l in open(f"{loc}den.dat",'r').readlines()][0][0])
  energy /= weight
  return { # Format matches result from read_afqmc for full AFQMC run.
      'warmup': 0,
      'safe_energy': energy,
      'energy': energy,
      'stdev': 0.0,
      'error': 0.0,
      #'blockbeta': [],
      #'blockdata': []
    }

class Hamiltonian:
  ''' Class for containing the Hamiltonian parts needed for a calculation.'''
  def __init__(self,onebody=None,twobody=None,nelec=None,constant=0.0):
    self.onebody = onebody
    self.twobody = twobody
    self.nelec = nelec
    self.constant = constant

  def to_hdf(self,hdf):
    ''' Export data as HDF5 file.'''
    with File(hdf,'w') as outf:
      outf.create_dataset("onebody",data=self.onebody)
      outf.create_dataset("twobody",data=self.twobody)
      outf.create_dataset("constant",data=self.constant)
      outf.create_dataset("nelec",data=self.nelec)

  def from_hdf(self,hdf):
    ''' Import data from HDF5 file.'''
    with File(hdf,'r') as inpf:
      self.onebody = inpf['onebody'][()]
      self.twobody = inpf['twobody'][()]
      self.nelec   = inpf['nelec'][()]
      self.constant = inpf['constant'][()]
    
    return self

if __name__=='__main__':
  main()
