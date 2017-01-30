import numpy as np
from inspect import getargspec
from scipy.optimize import curve_fit
import os
import matplotlib.pyplot as plt

# Nicer colors for plotting from colorbrewer2.org.
# Common colors are first letter only, prefix with l or d means light or dark.
pc = {
    # 12-Class paired.
    'lb':'#a6cee3', # Light blue.
    'b' :'#1f78b4', # Blue.
    'lg':'#b2df8a', # Light green.
    'g' :'#33a02c', # Green.
    'lr':'#fb9a99', # Light red.
    'r' :'#e31a1c', # Red
    'lo':'#fdbf6f', # Light orange/peach.
    'o' :'#ff7f00', # Orange.
    'lp':'#cab2d6', # Light purple/fushia.
    'p' :'#6a3d9a', # Purple.
    'lt':'#ffff99', # Light tan
    'tan' :'#b15928', # Cowhide tan.
    # 8-Class dark.
    't':    '#1b9e77', # Teal.
    'do':   '#d95f02', # Dark orange.
    'pink': '#e7298a', # Dark pink.
    'dy':   '#e6ab02', # Dark yellow/gold.
    'dgray':'#666666', # Dark gray.
    # 9-class Set1.
    'dr':    '#e41a1c', # Dark red.
    'db':    '#377eb8', # Dark blue.
    'y':     '#ffff33', # Yellow.
    'brown': '#a65628', 
    'gray':  '#999999'
  }

# Sets of colors to automatically choose from.
ps = {
    'dark8':['#0b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666'],
    'cb12':['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']
  }
pm = ["o",
    "s",
    "p",
    "h",
    "H",
    "D",
    "d"
    "v",
    "^",
    "<",
    ">",
    "*",
    ".",
    ",",
    "1",
    "2",
    "3",
    "4",
    "8",
    "+",
    "x",
]

# My own plotting defaults.
myplotdef={
    'mew':0.5,
    'mec':'k'
  }
myerrdef={
    'capthick':1,
    'capwidth':1
  }

notes = """
Things I commonly have to look up:
  For reference: all the marker choices for matplotlib:
   "."         point
   ","         pixel
   "o"         circle
   "v"         triangle_down
   "^"         triangle_up
   "<"         triangle_left
   ">"         triangle_right
   "1"         tri_down
   "2"         tri_up
   "3"         tri_left
   "4"         tri_right
   "8"         octagon
   "s"         square
   "p"         pentagon
   "*"         star
   "h"         hexagon1
   "H"         hexagon2
   "+"         plus
   "x"         x
   "D"         diamond
   "d"         thin_diamond
   "|"         vline
   "_"         hline
   TICKLEFT    tickleft
   TICKRIGHT   tickright
   TICKUP      tickup
   TICKDOWN    tickdown
   CARETLEFT   caretleft
   CARETRIGHT  caretright
   CARETUP     caretup
   CARETDOWN   caretdown
   "None"      nothing
   None        nothing
   " "         nothing
   ""          nothing
Color maps:
  Good ``uniform'' sequential: viridis,plasma.
  Good diverging: Spectral,seismic,bwr,BrBG
"""

def matplotlib_header(usetex=True,family='serif'):
  import seaborn as sns
  sns.set_style("ticks")
  plt.rc('axes.formatter',useoffset=False)
  plt.rc('text',usetex=usetex)
  plt.rc('font',style='normal')
  plt.rc('font',family=family)
  plt.rc('font',serif='Computer Modern')
  plt.rc('text',usetex=True)
  ticksize = 4
  plt.rc('xtick.major',size=ticksize)
  plt.rc('ytick.major',size=ticksize)

def fix_lims(ax_inp,factor=0.04):
  """
  Create buffer around all points that is factor*(data range) wide.
  """
  minx,miny = np.Inf,np.Inf
  maxx,maxy = -np.Inf,-np.Inf
  if type(ax_inp) != np.ndarray:
    ax_array = np.array((ax_inp))
  else:
    ax_array = ax_inp
  for ax in ax_array.flatten():
    for line in ax.get_lines():
      # axvlines store the data as lists and often should be ignored.
      if type(line.get_xdata()) is type([]):
        continue
      if line.get_xdata().max() > maxx:
        maxx = line.get_xdata().max()
      if line.get_ydata().max() > maxy:
        maxy = line.get_ydata().max()
      if line.get_xdata().min() < minx:
        minx = line.get_xdata().min()
      if line.get_ydata().min() < miny:
        miny = line.get_ydata().min()
  xs = factor*(maxx-minx)
  ys = factor*(maxy-miny)
  for ax in ax_array.flatten():
    if xs != 0: ax.set_xlim(minx-xs,maxx+xs)
    if ys != 0: ax.set_ylim(miny-ys,maxy+ys)

def slope(x,y): return (y[-1]-y[0])/(x[-1]-x[0])

def thin_ticks(ticks,div=2,start=0,shift=0,append=0):
  newticks = [ticks[div*i-shift] for i in range(start,len(ticks)//div+append)]
  return newticks

def idxmap(arraylike):
  """ Map elements of an array to it's index (reverse mapping of array itself)."""
  return dict(zip(arraylike,range(len(arraylike))))

class FitFunc:
  """
  Define a function that can be fit to and plotted with.
  """
  def __init__(self,form,jacobian=None,pnames=[]):
    self.form = form
    self.jac  = jacobian
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

  def fit(self,xvals,yvals,evals=None,guess=(),handle_nans=True,**kwargs):
    """
    Use xvals and yvals +/- evals to fit params with initial values p0.

    evals == None means don't use errorbars.
    guess == () means guess all 1.0 for the parameters (usually bad!)
    kwargs passed to curve_fit()
    handle_nans automatically drops tuples that have nan in any of xvals, yvals,
      or evals.
    """
    if handle_nans:
      drop = np.isnan(xvals)
      drop = drop | np.isnan(yvals)
      if evals is not None:
        drop = drop | np.isnan(evals)
      xvals = np.array(xvals)[~drop]
      yvals = np.array(yvals)[~drop]
      if evals is not None:
        evals = np.array(evals)[~drop]
    if guess == ():
      guess = self._set_default_parms(xvals,yvals,evals)
    if (evals is None) or (np.isnan(evals).any()):
      fit = curve_fit(self.form,
        xvals,yvals,
        p0=guess,**kwargs)
    else:
      fit = curve_fit(self.form,
        xvals,yvals,sigma=evals,
        absolute_sigma=True,
        p0=guess,**kwargs)
    self.parm = np.array(guess)
    self.perr = np.array(guess)
    for pi,p in enumerate(guess):
      self.parm[pi] = fit[0][pi]
      self.perr[pi] = fit[1][pi][pi]**.5
    self.cov  = fit[1]
    if len(self.parm) == len(self.pnms):
      self.pmap = dict(zip(self.pnms,self.parm))
      self.emap = dict(zip(self.pnms,self.perr))

  def eval(self,x):
    """
    Evaluate fitted function at point x.
    """
    if self.parm is None: return None
    else:
      return self.form(x,*self.parm)

  def eval_error(self,x):
    """
    Error from evaluating fitted function at point x.
    """
    if (self.parm is None) or (self.perr is None) or (self.jac is None): return None
    else:
      return np.dot( self.jac(x,*self.parm).T,
                     np.dot(self.cov,
                            self.jac(x,*self.parm)))**.5

  def get_parm(self,key="print"):
    if key=="print":
      print("What did you want? Available keys:")
      return ', '.join(self.pmap.keys())
    elif self.pmap != {}:
      return self.pmap[key]
    else:
      print("You must set pnames to use get_parm().")
      print("Alternatively, use self.parm.")
      return None

  def get_parm_err(self,key="print"):
    if key=="print":
      print("What did you want? Available keys:")
      return ', '.join(self.pmap.keys())
    elif self.pmap != {}:
      return self.emap[key]
    else:
      print("You must set pnames to use get_parm().")
      print("Alternatively, use self.perr.")
      return None

  def _set_default_parms(self,xvals,yvals,evals):
    # Better things possible for more specific functions.
    return np.ones(len(getargspec(self.form)[0])-1)

class LinearFit(FitFunc):
  """
  FitFunc of form c*x + y0
  """
  def __init__(self,pnames=['slope','yint']):
    def form(x,c,y0):
      return c*x + y0
    def jac(x,c,y0):
      return np.array([x,1.0]).T
    self.form = form
    self.jac  = jac
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

  def _set_default_parms(self,xvals,yvals,evals):
    return (slope(xvals,yvals), yvals[abs(xvals).argmin()])

class LinearFit_xcross(FitFunc):
  """
  FitFunc of form c*(x - x0)
  """
  def __init__(self,pnames=['slope','xint']):
    def form(x,c,x0):
      return c*(x - x0)
    def jac(x,c,y0):
      return np.array([x,-c]).T
    self.form = form
    self.jac  = jac
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

  def _set_default_parms(self,xvals,yvals,evals):
    return (slope(xvals,yvals), xvals[abs(yvals).argmin()])

class QuadraticFit(FitFunc):
  """
  FitFunc of form c*(x - xm)**2 + yc
  """
  def __init__(self,pnames=['quadratic','xmin','ycrit']):
    def form(x,c,xm,yc):
      return c*(x - xm)**2 + yc
    def jac(x,c,xm,yc):
      return np.array([(x-xm)**2,2*c*(x-xm),1]).T
    self.form = form
    self.jac  = jac
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

class CubicFit(FitFunc):
  """
  FitFunc of form a*(x - xm)**3 + b*(x - xm)**2 + yc
  """
  def __init__(self,pnames=['cubic','quadratic','xmin','ycrit']):
    def form(x,a,b,xm,yc):
      return a*(x - xm)**3 + b*(x - xm)**2 + yc
    def jac(x,a,b,xm,yc):
      return np.array([(x-xm)**3,(x-xm)**2,-3*a*(x-xm)**2-2*b*(x-xm),1]).T
    self.form = form
    self.jac  = jac
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

  def _set_default_parms(self,xvals,yvals,evals):
    # These will work well for cubics that are close to parabolic, with samples
    # centered around the min or max.
    if yvals[yvals.shape[0]//2] > yvals[0]:
      return (1.0, -1.0, xvals.mean(), yvals.max())
    elif yvals[yvals.shape[0]//2] < yvals[0]:
      return (1.0, 1.0, xvals.mean(), yvals.min())
    else: # It's something flat-ish?
      return (0.0, 0.0, xvals.mean(), yvals.min())

class CubicFit_fixmin(FitFunc):
  """
  FitFunc of form a*(x - xm)**3 + b*(x - xm)**2 + yc
  """
  def __init__(self,xm,pnames=['cubic','quadratic','ycrit']):
    def form(x,a,b,yc):
      return a*(x - xm)**3 + b*(x - xm)**2 + yc
    def jac(x,a,b,yc):
      return None # implement me!
    self.form = form
    self.jac  = jac
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

class CubicFit_zeros(FitFunc):
  """
  FitFunc of form a(x-x1)(x-x2)(x-x3)
  """
  def __init__(self,pnames=['cubic','zero1','zero2','zero3']):
    def form(x,a,x1,x2,x3):
      return a*(x-x1)*(x-x2)*(x-x3)
    def jac(x,a,b,yc):
      return None # implement me!
    self.form = form
    self.jac  = jac
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

class NormedGaussianFit(FitFunc):
  """
  FitFunc of form (2 pi sigma)**-0.5 exp(-(x-mu)**2/2sigma**2)
  """
  def __init__(self,pnames=['mean','std']):
    def form(x,mu,sigma):
      return (2.*np.pi*sigma)**-0.5 * np.exp(-(x-mu)**2/2./sigma**2)
    def jac(x,mu,sigma):
      return None 
    self.form = form
    self.jac  = None # implement me!
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

  def _set_default_parms(self,xvals,yvals,evals):
    return (np.mean(xvals),np.std(xvals))

class GaussianFit(FitFunc):
  """
  FitFunc of form A/(2 pi sigma)**0.5 exp(-(x-mu)**2/2sigma**2)
  """
  def __init__(self,pnames=['mean','std','amp']):
    def form(x,mu,sigma,A):
      return A/(2.*np.pi*sigma)**0.5 * np.exp(-(x-mu)**2/2./sigma**2)
    def jac(x,mu,sigma):
      return None 
    self.form = form
    self.jac  = None # implement me!
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

  def _set_default_parms(self,xvals,yvals,evals):
    return (np.mean(xvals),np.std(xvals),max(yvals)-min(yvals))

class EOSFit(FitFunc):
  """
  Anton-Shmidt (DOI: 10.1016/S0966-9795(97)00017-4) Equation of state E(V):
  P(V) = -b(V/V0)^n log(V/V0)
  => E(V) = bV0/(n+1) (V/V0)^(n+1) (ln(V/V0) - 1/(n+1)) + Einf
  """
  def __init__(self,pnames=['bulk_mod','eq_vol','n','Einf']):
    def energy(V,b,V0,n,Einf):
      return b*V0/(n+1) * (V/V0)**(n+1) * (np.log(V/V0) - 1/(n+1)) + Einf
    def pressure(V,b,V0,n):
      return -b*(V/V0)**n  * np.log(V/V0)

    self.form = energy
    self.derv = pressure
    self.jac  = None # Haven't bothered yet.
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

  def eval_derv(self,x):
    """
    Evaluate derivative function at point x.
    """
    if self.parm is None: return None
    else:
      return self.derv(x,*self.parm[:-1])

  def _set_default_parms(self,xvals,yvals,evals):
    # Should be good if you're only after positive pressures.
    HaA3_GPa = 4359.74434
    maxyidx = yvals.argmax()
    return (
        1.0/HaA3_GPa,   # b: rough scale for bulk modulus.
        xvals[maxyidx], # V0: Volume at ambient pressure.
        -2.0,           # n: suggested by Anton et al.
        yvals[maxyidx]  # Einf: rough energy scale near V0.
      )

class EOSFit_fixV0(EOSFit):
  def __init__(self,V0,pnames=['bulk_mod','n','Einf']):
    def energy(V,b,n,Einf):
      return b*V0/(n+1) * (V/V0)**(n+1) * (np.log(V/V0) - 1/(n+1)) + Einf
    def pressure(V,b,n):
      return -b*(V/V0)**n  * np.log(V/V0)

    self.form = energy
    self.derv = pressure
    self.jac  = None # Haven't bothered yet.
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

class EOSFit_fixn(EOSFit):
  def __init__(self,n,pnames=['bulk_mod','V0','Einf']):
    def energy(V,b,V0,Einf):
      return b*V0/(n+1) * (V/V0)**(n+1) * (np.log(V/V0) - 1/(n+1)) + Einf
    def pressure(V,b,V0):
      return -b*(V/V0)**n  * np.log(V/V0)

    self.form = energy
    self.derv = pressure
    self.jac  = None # Haven't bothered yet.
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

class EOSFit_fixV0_fixn(EOSFit):
  def __init__(self,V0,n,pnames=['bulk_mod','Einf']):
    def energy(V,b,Einf):
      return b*V0/(n+1) * (V/V0)**(n+1) * (np.log(V/V0) - 1/(n+1)) + Einf
    def pressure(V,b):
      return -b*(V/V0)**n  * np.log(V/V0)

    self.form = energy
    self.derv = pressure
    self.jac  = None # Haven't bothered yet.
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

class MorseFit(FitFunc):
  """
  Morse potential for model of bonding.
  V(r) = Deq(1 - exp(-a(r-req)))^2 + V(req)

  Suggested guesses: 
  Deq  = max(V) - min(V)
  req  = argmin(V)
  a    ~ 0.05
  Vreq = min(V)
  """
  def __init__(self,pnames=['depth','exp_coef','eq_radius','eq_potential']):
    def pot(r,Deq,a,req,Vreq):
      return Deq*(1. - np.exp(-a*(r-req)))**2 + Vreq
    self.form = pot
    self.jac  = None # Haven't bothered yet.
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None
  def _set_default_parms(self,xvals,yvals,evals):
    # Should be good if you're only after positive pressures.
    minr = yvals.argmin()
    return (
        yvals[-1] - yvals[minr],# Depth = Asytote - min.
        xvals[minr],            # Position of minimum.
        0.05,                   # Emperical suggestion.
        yvals[minr]             # Bottom of potential.
      )

# Its been a while and I'm not sure how this is different. Delete?
class MorseFitpp(FitFunc):
  """
  Morse potential for model of bonding.
  V(r) = Deq(1 - exp(-a(r-req)))^2 + V(req)

  Suggested guesses: 
  Deq  = max(V) - min(V)
  req  = argmin(V)
  a    ~ 0.05
  Vreq = min(V)
  """
  def __init__(self,pnames=['depth','exp_coef','eq_radius','eq_potential']):
    def pot(r,Deq,a,req,Vreq):
      return Deq*(1. - np.exp(-a*(r-req)**3))**2 + Vreq
    self.form = pot
    self.jac  = None # Haven't bothered yet.
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None

class CatagoryPlot: 
  """ Use a pandas DataFrame to make plots broken down by color, row, column,
  and marker. Somewhat similar to what ggplot can handle (more elegantly).

  For example: df.columns=[x,y,z],
  cp=CatagoryPlot(df,color='z',mark='z')
  cp.plot('x','y')

  Now cp.fig will have a figure of x vs y with color and marker broken down by z.
  
  Call self.plot() to actually make a plot. 
  
  There are a lot of bugs I know about but haven't bothered to fix. If you find
  one, I can probably fix it pretty quick, so long as there's a desire to have
  it fixed.
  """
  def __init__(self,df,
      row='dummy',col='dummy',
      color='dummy',mark='dummy',
      cmap=None,mmap=None,sharex=False,sharey=False):

    if 'dummy' in df.columns:
      print("CatagoryPlot: Warning, I'm not going to use the 'dummy' column!")


    assert df.shape[0]>0 and df.shape[1]>0, "Empty dataframe!"
    self.fulldf=df
    self.fulldf['dummy']='dummy'
    self.row=row
    self.col=col
    self.color=color
    self.mark=mark
    self.plotargs={}
    if cmap is None:
      unique_cols=self.fulldf[color].unique()
      self.cmap=dict(zip(unique_cols,ps['dark8'][:unique_cols.shape[0]]))
    else: 
      self.cmap=cmap
    if mmap is None:
      unique_marks=self.fulldf[mark].unique()
      self.mmap=dict(zip(unique_marks,pm[:unique_marks.shape[0]]))
    else: 
      self.cmap=cmap

    self.fig,self.axes=plt.subplots(
        self.fulldf[row].unique().shape[0],
        self.fulldf[col].unique().shape[0],
        squeeze=False,sharex=sharex,sharey=sharey
      )
    self.rowmap=idxmap(df[row].unique())
    self.colmap=idxmap(df[col].unique())

  def plot(self,xvar,yvar,plotargs={}):
    self.plotargs=plotargs
    for lab,df in self.fulldf.groupby([self.row,self.col,self.mark,self.color]):
      self.axes[self.rowmap[lab[0]],self.colmap[lab[1]]]\
          .plot(df[xvar],df[yvar],self.mmap[lab[2]],
          color=self.cmap[lab[3]],**plotargs)
      self.fig.tight_layout()

  def add_legend(self,axidx=(0,0),labstr="%s,%s",labmap={},args={}):
    """ Make a legend for the markers and/or colors. labstr controls how the
    labels combine the two parameters (if they're distict). labmap maps data to
    pretty labels. locargs is passed to axes.legend(). """
    unique_cols=self.fulldf[self.color].unique()
    unique_marks=self.fulldf[self.mark].unique()
    if self.mark=='dummy': 
      if labmap=={}:
        labmap=dict(zip(unique_cols,unique_cols))
      prox=[plt.Line2D([],[],
            linestyle='',
            marker=self.mmap['dummy'],color=self.cmap[unique],label=labmap[unique],
            **self.plotargs
          ) for unique in unique_cols
        ]
    elif self.color=='dummy': 
      if labmap=={}:
        labmap=dict(zip(unique_marks,unique_marks))
      prox=[plt.Line2D([],[],
            linestyle='',
            marker=self.mmap[unique],color=self.cmap['dummy'],label=labmap[unique],
            **self.plotargs
          ) for unique in unique_cols
        ]
    else:
      if labmap=={}:
        labmap=dict(zip(unique_marks,unique_marks))
        labmap.update(dict(zip(unique_cols,unique_cols)))
      prox=[plt.Line2D([],[],
            linestyle='',
            marker=self.mmap[unique_mark],
            color=cmap[unique_col],
            label=labstr%(labmap[unique_col],labmap[unique_mark]),
            **self.plotargs
          )
          for unique_col in unique_cols for unique_mark in unique_marks
        ]
    self.axes[axidx].legend(handles=prox,loc='best')
