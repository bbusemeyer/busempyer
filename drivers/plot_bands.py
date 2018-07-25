import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plot_tools as pt
from read_bands import read_fort25
import pickle as pkl
ev=27.2114
pt.matplotlib_header()

def plot_bands(spin=1):
  dat = read_fort25("fort.25")
  totk=0.0
  nbands=dat['bands'][0]['dat'].shape[1]
  print(nbands)
  fullbands=[[] for i in range(nbands)]
  fullkpts=[]
  breaks=[0.0]
  klabs=['000']
  for leg in dat['bands'][:len(dat['bands'])//1]:
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
  for bidx,band in enumerate(fullbands):
    if -3<band[0]<3:
      print(bidx)
      c='r'
      l=2
    else:
      c='k'
      l=1
    ax.plot(fullkpts,band,c,lw=l)
  ax.set_xticks(breaks)
  ax.set_xticklabels(klabs)
  ax.set_xlim((min(breaks)-0.1,max(breaks)+0.1))
  #ax.set_ylim((-4,4))
  for line in breaks:
    ax.axvline(line,color='k',lw=0.5)

  #models=pkl.load(open('../../analysis/kptdf.pkl','rb'))
  #moddf=models.apply(lambda x:x.params).drop('020').reset_index()
  #print(moddf)
  #moddf['kplace']=moddf['kpoint'].apply(lambda x: breaks[klabs.index(x.replace('2','1'))])
  #bandcols=["e_{%s,%s}"%(band,band) for band in range(8)]
  #shift=dat['bands'][0]['dat'][0,14]-dat['bands'][0]['efermi']-moddf.loc[moddf['kpoint']=='000','e_{4,4}'].squeeze()
  #moddf[bandcols]=(moddf[bandcols]+shift)*ev
  #extra=moddf[moddf['kpoint']=='000'].copy()
  #extra['kplace']=breaks[-1]
  #moddf=pd.concat([moddf,extra])

  #modelbands*=ev

  #for band in bandcols:
  #  ax.plot(moddf['kplace'],moddf[band],'d',color=pt.pc['b'])
 
  ax.set_xlabel(r"$k$-point$(\frac{k_x}{\pi a}\frac{k_y}{\pi b}\frac{k_z}{\pi c})$")
  ax.set_ylabel("$E - E_\mathrm{F}$ (eV)")
  fig.set_size_inches(6,3)
  fig.tight_layout()
  fig.savefig("bands.pdf")
  fig.savefig("bands.eps")

if __name__=='__main__':
  plot_bands()
