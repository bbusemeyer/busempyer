import numpy as np

def read_band_section(inp):
  line = inp.readline().split()
  nkpt = int(line[line.index('nkpt=') + 1])
  kpts, evals = [],[]
  for kpt in range(nkpt):
    line  = inp.readline().split()
    nband = int(line[line.index('nband=')+1][:-1])
    p     = line.index('kpt=') + 1
    kpts.append( np.array(line[p:p+3],dtype=float) )

    newval = []
    for bi in range((nband-1)/8 + 1): 
      newval += inp.readline().split() 
    evals.append( np.array(newval,dtype=float) )
    
  plen = [((kpts[i]-kpts[i-1])**2).sum() for i in range(1,len(kpts))]
  plen = np.array([0.0] + plen)
  plen = plen.cumsum()
  return [np.array(kpts),plen,np.array(evals)]
