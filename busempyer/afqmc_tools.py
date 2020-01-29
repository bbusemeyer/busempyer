''' Misc tools for interfacing with AFQMCLab.'''
import pyblock, numpy, os
from pandas import DataFrame

def main():
  print("No default actions.")

def dump_afqmc_params(**opts):
  ''' dump opts to afqmc_params for AFQMCLab.'''
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
    return {}

  edf = DataFrame([l.split() for l in open(f"{loc}HNum.dat",'r').readlines()],
      columns=('energy','imenergy'),dtype=float)
  edf = edf.join(DataFrame([l.split() for l in open(f"{loc}den.dat",'r').readlines()],
      columns=('weight','imweight'),dtype=float))
  edf = edf.join(DataFrame([l.split() for l in open(f"{loc}beta.dat",'r').readlines()],
      columns=['beta'],dtype=float))
  edf['beta'] = edf['beta'].round(10)
  assert edf['beta'].shape == edf['beta'].unique().shape, "Multiple betas: did you make a dirty restart?"
  safe_energy = edf['energy'].sum() / edf['weight'].sum()
  edf['energy'] = edf['energy'] / edf['weight']

  if warmup is None:
    start,count = 0,0
    for start in range(edf.shape[0]):
      if edf.loc[start,'energy'] < safe_energy:
        count += 1
        if count >= 5: break
    warmup = start

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

if __name__=='__main__':
  main()
