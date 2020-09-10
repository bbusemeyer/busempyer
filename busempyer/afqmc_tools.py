''' Misc tools for interfacing with AFQMCLab.'''
import pyblock, numpy, os
from pandas import DataFrame
from busempyer.qmcdata import estimate_warmup
from h5py import File

def main():
  print("No default actions.")

def dump_afqmc_param(**opts):
  ''' dump opts to afqmc_param for AFQMCLab.'''
  outlines = [f"{key} {opts[key]}" for key in opts]
  with open('afqmc_param','w') as outf:
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
        'blockbeta': None,
        'blockdata': None
      }

  blockenergy = pyblock.blocking.reblock(edf.iloc[warmup:]['energy'].values)
  optimal = pyblock.blocking.find_optimal_block(edf.shape[0]-warmup,blockenergy)[0]
  get = optimal if not numpy.isnan(optimal) else 0
  if get==0: print("Warning: no optimal block found, errorbars will be inaccurate.")
  ndata = blockenergy[get].ndata
  earr = edf.iloc[warmup:]['energy'].values
  earr = earr[earr.shape[0]%ndata:].reshape(ndata,earr.shape[0]//ndata)
  blockdf = DataFrame({
      'energy':earr.mean(axis=1),
      'stdev':earr.std(axis=1),
    })

  blocknbeta = (edf.shape[0]-warmup)//blockdf.shape[0]
  results =  {
      'warmup':       warmup,
      'safe_energy':  safe_energy,
      'energy':       blockdf['energy'].mean(),
      'stdev':        blockdf['energy'].std(),
      'error':        blockdf['energy'].std()/blockdf.shape[0]**0.5,
      'blockbeta':    (edf.iloc[warmup:]['beta'].values[blocknbeta//2::blocknbeta]).tolist(),
      'blockdata':    blockdf['energy'].values.tolist(),
    }

  if return_trace:
    results['beta']  = edf['beta'].values.tolist()
    results['trace'] = edf['energy'].values.tolist()

  return results 

def read_raw_afqmc(loc="./"):
  edf = DataFrame([l.split() for l in open(f"{loc}HNum.dat",'r').readlines()],
      columns=('energy','imenergy'),dtype=float)
  edf = edf.join(DataFrame([l.split() for l in open(f"{loc}den.dat",'r').readlines()],
      columns=('weight','imweight'),dtype=float))
  edf = edf.join(DataFrame([l.split() for l in open(f"{loc}beta.dat",'r').readlines()],
      columns=['beta'],dtype=float))
  edf['beta'] = edf['beta'].round(10)
  assert edf['beta'].shape == edf['beta'].unique().shape, "Multiple betas: did you make a dirty restart?"

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
      'warmup': None,
      'safe_energy': energy,
      'energy': energy,
      'stdev': 0.0,
      'error': 0.0,
      'blockbeta': [],
      'blockdata': []
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
