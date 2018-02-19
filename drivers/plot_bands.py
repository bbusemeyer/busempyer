import matplotlib.pyplot as plt
import numpy as np
import plot_tools as pt
from read_bands import read_fort25
ev=27.2114
pt.matplotlib_header()

def plot_bands():
  dat = read_fort25("fort.25")
  totk=0.0
  nbands=dat['bands'][0]['dat'].shape[1]
  fullbands=[[] for i in range(nbands)]
  fullkpts=[]
  breaks=[0.0]
  klabs=['000']
  for leg in dat['bands'][:len(dat['bands'])//2]:
    bands=(leg['dat']-leg['efermi'])*ev
    dk=leg['dkp']

    kpts=totk+dk*np.arange(bands.shape[0])
    fullkpts+=kpts.tolist()
    breaks.append(kpts[-1])
    klabs.append(leg['k1'])
    for bi,band in enumerate(bands.T):
      fullbands[bi]+=band.tolist()
    totk+=dk*bands.shape[0]

  fig,ax=plt.subplots(1,1)
  for band in fullbands:
    if -3<band[0]<3:
      c='r'
      l=2
    else:
      c='k'
      l=1
    ax.plot(fullkpts,band,c,lw=l)
  ax.set_xticks(breaks)
  ax.set_xticklabels(klabs)
  ax.set_xlim((min(breaks),max(breaks)))
  ax.set_ylim((-4,4))
  for line in breaks:
    ax.axvline(line,color='k',lw=0.5)

  
  ax.set_xlabel(r"$k$-point$(\frac{k_x}{\pi a}\frac{k_y}{\pi b}\frac{k_z}{\pi c})$")
  ax.set_ylabel("$E - E_\mathrm{F}$ (eV)")
  fig.set_size_inches(3,3)
  fig.tight_layout()
  fig.savefig("bands.pdf")
  fig.savefig("bands.eps")

if __name__=='__main__':
  plot_bands()
