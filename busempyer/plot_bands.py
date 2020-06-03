''' A code fro plotting bands from CRYSTAL.'''
from matplotlib.pyplot import subplots
from numpy import arange
from busempyer.read_bands import read_fort25
from busempyer.plot_tools import matplotlib_header
matplotlib_header()

EV = 27.2114
KMAP = {
    '000':'$\Gamma$',
    '100':'X',
    '110':'M',
    '101':'R',
    '111':'A',
    '001':'Z'
  }

def main():
  # Sample useage:
  #erange = (-5,3)
  #plot_bands('unpolarized/fort.25',erange,spin=0,figname='../plots/fu4/bands_unp')
  #plot_bands('checkerboard_broy/fort.25',erange,spin=1,figname='../plots/fu4/bands_che')
  #plot_bands('collinear/fort.25',erange,spin=1,figname='../plots/fu4/bands_col')
  print('No default usage')

def plot_bands(loc,lims,spin,red=(0,0),figname=None):
  ''' Plots bands from loc.
  Args:
    loc (str): path to fort.25 file.
    lims (tuples): Energy limits for plot.
    spin (int): 0=all (unpolarized), 1=first half (spin up), 2=second half (spin down).
    red (tuple): range of bands to color red.
    figname (str): root name for file to export plot into.
  Returns:
    fig,ax: Figure and axes object with plot.
  '''
  fort25 = read_fort25(loc)

  #print(fort25['bands'][0].keys())
  #print(fort25['bands'][0]['dat'].shape)
  #print(fort25['bands'][0]['k0'])
  #print(fort25['bands'][0]['k0'])
  #print(fort24['bands'][-1]['dkp'])

  fig,ax = subplots(1,1)

  if spin==1:
    want_bands = fort25['bands'][:len(fort25['bands'])//2] 
  elif spin==2:
    want_bands = fort25['bands'][len(fort25['bands'])//2:] 
  else:
    want_bands = fort25['bands']

  starts = [0]
  for leg in want_bands:
    nkp = leg['dat'].shape[0]
    dktot = nkp*leg['dkp']
    for bidx,band in enumerate(leg['dat'].T):
      if bidx in range(*red): color='r'
      else: color = 'k'
      ax.plot(starts[-1] + leg['dkp']*arange(nkp),(band - leg['efermi'])*EV,color+'-',lw=1,alpha=0.5)
    starts.append(starts[-1]+dktot)

  for start in starts[1:-1]:
    ax.axvline(start,color='k',lw=1)
  ax.axhline(0.0,color='k',lw=1)
  ax.set_xlim(starts[0],starts[-1])
  ax.set_ylim(lims)
  ax.set_xticks(starts)
  ax.set_xticklabels([KMAP[leg['k0']] for leg in want_bands] + [KMAP[want_bands[-1]['k1']]])
  ax.set_xlabel('Kpoint')
  ax.set_ylabel('Bands (eV)')
  fig.set_size_inches(4,3)
  fig.tight_layout()

  if figname is not None:
    print("Saving %s."%figname)
    fig.savefig(figname+'.pdf')
    fig.savefig(figname+'.eps')
    fig.savefig(figname+'.png',dpi=400)
  
  return fig,ax

def compare_pyscf():
  import sys
  sys.path.append('..')
  from make_iaos import load_pyscf

  fort25 = read_fort25('collinear/fort.25')
  bands = fort25['bands'][0]['dat'][0,:]

  cell,mf = load_pyscf('collinear/GRED.DAT','collinear/KRED.DAT','collinear/crys.o')
  gamma = mf.mo_energy[0,0,:bands.shape[0]]

  print( abs(gamma - bands).max() )

if __name__=='__main__':
  main()
