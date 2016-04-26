import numpy as np
from inspect import getargspec
from scipy.optimize import curve_fit
import os

def fix_lims(ax_array,factor=0.04):
  """
  Create buffer around all points that is factor*(data range) wide.
  """
  minx,miny = np.Inf,np.Inf
  maxx,maxy = -np.Inf,-np.Inf
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

  def fit(self,xvals,yvals,evals=None,guess=(),**kwargs):
    """
    Use xvals and yvals +/- evals to fit params with initial values p0.

    evals == None means don't use errorbars.
    guess == () means guess all 1.0 for the parameters (usually bad!)
    kwargs passed to curve_fit()
    """
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
