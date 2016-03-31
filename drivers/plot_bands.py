import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
ticksize = 4
sns.set_style("ticks")
plt.rc('xtick.major',size=ticksize)
plt.rc('ytick.major',size=ticksize)
plt.rc('font',style='normal')
plt.rc('font',family='serif')
plt.rc('font',serif='Computer Modern Roman')
plt.rc('text',usetex=True)
plt.rc('lines',linewidth=1)
plt.rc('axes.formatter',useoffset=False)
from read_bands import read_fort25
ev = 27.211396132

dat = read_fort25("fort.25")
bands = (dat['bands'][0]['dat'] - dat['bands'][0]['efermi'])*ev*1e3
dk    = dat['bands'][0]['dkp']

fig,ax = plt.subplots()
fig.set_size_inches(3,3)
for band in bands.T:
  full = band.tolist()
  kpts = [dk*(i - len(full)) for i in range(2*len(full))]
  full = list(reversed(full))+full
  ax.plot(kpts,full,color='k')
ax.set_ylabel("$E - E_\mathrm{F}$ (meV)")
#ax.set_xlim(-0.2,0.2)
#ax.set_xticks(np.linspace(-0.1,0.1,3))
ax.set_ylim(-500,800)
fig.tight_layout()
fig.savefig("unp_bands.pdf")
