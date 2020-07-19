#!/usr/bin/env python3
''' Analysis for QMC trace.'''
from numpy import isnan
from pyblock.blocking import find_optimal_block,reblock

def main():
  print("No main functionality.")

def estimate_warmup(energy):
  ''' Estimate equillibration of data assuming that at least half the data is equillibrated.'''
  safe = energy[energy.shape[0]//2:]

  safeblocks = reblock(safe)
  bo = find_optimal_block(safe.shape[0],safeblocks)[0]
  if isnan(bo): 
    assert 0, "Not sure what should happen here!"
    bo = -1
  blockdata = safeblocks[bo]

  blocks = energy[energy.shape[0]%(2*blockdata.ndata):].reshape((2*blockdata.ndata),energy.shape[0]//(2*blockdata.ndata))
  blocks = blocks.mean(axis=1)

  signflags = [False,False]
  for bi,block in enumerate(blocks):
    signflags[0] |= block >= blockdata.mean
    signflags[1] |= block <= blockdata.mean
    if signflags[0] and signflags[1]: break

  return energy.shape[0]//(blockdata.ndata*2) * bi

if __name__=='__main__': main()
