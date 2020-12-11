#!/usr/bin/env python3
''' Analysis for QMC trace.'''
from pandas import DataFrame
from numpy import isnan
from pyblock.blocking import find_optimal_block,reblock

def main():
  print("No main functionality.")

def estimate_warmup(energy):
  ''' Estimate equillibration of data assuming that at least half the data is equillibrated.'''
  safe = energy[energy.shape[0]//2:]

  safeblocks = reblock(safe)
  try:
    bo = find_optimal_block(safe.shape[0],safeblocks)[0]
  except IndexError:
    print("\n !!! Error in blocking: maybe there is too little statistics. !!! ")
    return -1
  if isnan(bo): 
    print("\n !!! Insufficient data for warm-up estimation! I will return -1 as a flag for this error. !!!")
    return -1
  blockdata = safeblocks[bo]

  blocks = energy[energy.shape[0]%(2*blockdata.ndata):].reshape((2*blockdata.ndata),energy.shape[0]//(2*blockdata.ndata))
  blocks = blocks.mean(axis=1)

  signflags = [False,False]
  for bi,block in enumerate(blocks):
    signflags[0] |= block >= blockdata.mean
    signflags[1] |= block <= blockdata.mean
    if signflags[0] and signflags[1]: break

  return energy.shape[0]//(blockdata.ndata*2) * bi

def block_data(data):
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
