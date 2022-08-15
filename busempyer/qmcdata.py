#!/usr/bin/env python3
''' Analysis for QMC trace.'''
from pandas import DataFrame
from math import ceil
from numpy import isnan
from qharv.reel.forlib.stats import corr as compute_autocorr

def main():
  print("No main functionality.")

def estimate_warmup(energy):
  ''' Estimate equillibration of data assuming that at least half the data is equillibrated.'''
  safe = energy[energy.shape[0]//2:]

  try:
    autotime = ceil(compute_autocorr(safe))
  except ValueError:
    print("Cannot compute warmup time, returning -1 as a flag.")
    return -1
  blocks = reblock(energy, autotime)
  safemean = safe.mean()

  signflags = [False,False]
  for bix,block in enumerate(blocks):
    if not signflags[0] and block >= safemean:  signflags[0] = True
    if not signflags[1] and block <= safemean:  signflags[1] = True
    if signflags[0] and signflags[1]:           break

  return autotime * bix

def reblock(data, blocklen):
  ''' Block-average data in chunks of blocklen (usually the autocorrelation time).'''
  stubpos = -(data.shape[0]%blocklen)
  if stubpos:
    data = data[:stubpos]
    stub = data[stubpos:]
  else:
    stub = None
  data = data.reshape(-1,blocklen).mean(axis=1)

  # What to do with stub? Seems a waste to throw it, why not average it with the final entry?
  if stub is not None:
    data[-1] = (data[-1] + stub.mean())/2

  return data

# Obsolete.
def block_data_pyblock(data):
  print("Warning: I've found this reblocking procedure to greatly overestimate the errors. Use block_data instead.")
  from pyblock.blocking import find_optimal_block,reblock
  ''' Takes serially-correlated data and returns independent sampled data.'''
  blocked = reblock(data)
  optimal = find_optimal_block(len(data),blocked)[0]
  get = optimal if not isnan(optimal) else 0
  if get==0: print("Warning: no optimal block found, errorbars will be inaccurate.")
  ndata = blocked[get].ndata
  earr = data[data.shape[0]%ndata:].reshape(ndata,data.shape[0]//ndata)
  blockdf = DataFrame({
      'value':earr.mean(axis=1),
      'stdev':earr.std(axis=1), # Might be useless.
    })

  # Note that Pandas uses correct std formula.

  return blockdf

if __name__=='__main__': main()
