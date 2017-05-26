def timeperstep_fromvmc(
    walltime,
    nblocks,
    steps_per_block,
    configs_per_proc=1):
  vmc_proc_time=walltime
  tot_vmc_configs=configs_per_proc*nblocks*steps_per_block
  timeperstep=vmc_proc_time/tot_vmc_configs
  print("Time per step:",timeperstep)
  return timeperstep

def required_time(std, target, timeperstep, timestep, nprocs=1):
  nblocks=(std/target)**2
  steps_per_block=1./timestep
  nsteps=steps_per_block*nblocks
  print('out',nsteps)
  time=timeperstep*nsteps
  print("Total time required: {} hours.".format(time/3600))
  walltime=time/nprocs 
  print("Total wall time given {} processors: {} hours.".format(nprocs,walltime/3600))
  return walltime 

if __name__=='__main__':
  # Lets test if this code is self consistent.
  import numpy as np
  time=np.random.randint(1000,10000)
  nprocs=np.random.randint(1,1000)
  nblocks=np.random.randint(1,100)
  steps_per_block=np.random.randint(1,100)
  nprocs=1
  print('inp',steps_per_block)

  std=np.random.rand()*10
  target=std/nblocks**0.5
  timestep=1./steps_per_block

  timeperstep=timeperstep_fromvmc(
      time,nblocks,steps_per_block,nprocs)

  required=required_time(std,target,timeperstep,timestep,nprocs)

  assert np.isclose(required,time),\
    "Functions are inconsistent! Times:\n {} {}".format(time,required)

  print("Test: functions are consistent!'\n {} vs {}".format(time,required))

