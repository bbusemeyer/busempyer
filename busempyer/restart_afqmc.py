#!/usr/bin/env python3
''' Set up a restart for an AFQMC job.'''
from busempyer.afqmc_tools import read_afqmc_param, dump_afqmc_param

def main():
  '''Command-line usage.'''
  from sys import argv
  path = argv[1] if len(argv)>1 else './'

  prep_restart(path)

def prep_restart(loc=".", updates=None):
  ''' Save the old afqmc_params and make a new one that restarts the job.'''
  params = read_afqmc_param(loc+"/afqmc_param")

  params['initialWalkerFlag'] = "readAllWalkers"
  params['seed'] = -1
  params['thermalSize'] = 0
  params['backGroundInit'] = "readFromFile"
  params['ETAdjustMaxSize'] = 0
  params['ETAndBackGroundGrowthEstimateMaxSize'] = 0
  if updates is not None:
    params.update(updates)
  
  with open(loc+"/past.afqmc_param",'a') as outf:
    with open(loc+"/afqmc_param",'r') as inpf:
      outf.write(inpf.read())

  dump_afqmc_param(fn=loc+'/afqmc_param', **params)

if __name__=='__main__': main()
