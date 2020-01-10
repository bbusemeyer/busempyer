''' Misc tools for interfacing with AFQMCLab.'''
import pandas, pyblock, numpy, os

def read_afqmc():
  if not os.path.exists('HNum.dat'):
    return {}

  edf = pandas.DataFrame([l.split() for l in open('HNum.dat','r').readlines()],
      columns=('energy','imenergy'),dtype=float)
  edf = edf.join(pandas.DataFrame([l.split() for l in open('den.dat','r').readlines()],
      columns=('weight','imweight'),dtype=float))
  edf = edf.join(pandas.DataFrame([l.split() for l in open('beta.dat','r').readlines()],
      columns=['beta'],dtype=float))
  edf['beta'] = edf['beta'].round(10)
  assert edf['beta'].shape == edf['beta'].unique().shape, "Multiple betas: did you make a dirty restart?"
  safe_energy = edf['energy'].sum() / edf['weight'].sum()
  edf['energy'] = edf['energy'] / edf['weight']

  blockenergy = pyblock.blocking.reblock(edf['energy'].values)
  optimal = pyblock.blocking.find_optimal_block(edf.shape[0],blockenergy)[0]
  get = optimal if not numpy.isnan(optimal) else 0
  if get==0: print("Warning: no optimal block found, errorbars will be inaccurate.")
  ndata = blockenergy[get].ndata
  earr = edf['energy'].values
  earr = earr[earr.shape[0]%ndata:].reshape(ndata,earr.shape[0]//ndata)
  blockdf = pandas.DataFrame({
      'energy':earr.mean(axis=1),
      'stdev':earr.std(axis=1),
    })

  return {
      'safe_energy':safe_energy,
      'energy':blockdf['energy'].mean(),
      'stdev':blockdf['energy'].std(),
      'error':blockdf['energy'].std()/blockdf.shape[0]**0.5,
      'trace':blockdf['energy'].values.tolist()
    }

def dump_afqmc_params(**opts):
  ''' format opts into string that can be written to afqmc_params.'''
  outlines = [f"{key} {opts[key]}" for key in opts]
  with open('afqmc_param','w') as outf:
    outf.write('\n'.join(outlines))
