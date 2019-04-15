''' Convenience tools for simple and quick fitting. '''

from inspect import getargspec
from scipy.optimize import curve_fit
import numpy as np

### Fitting tools.

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

class LogFit(FitFunc):
  ''' 
  Logarithmic fit of the form:

  y=y0 + b log(x0-x)

  Note that if x0=0 for you, then you should just use the linear fit!
  '''
  def __init__(self,pnames=['intercept','slope','shift']):
    def logform(x,y0,b,x0):
      return y0 + b*np.log(x0-x)
    def logjac(x,y0,b,x0):
      return np.array((1,np.log(x0-x),-b/(x0-x)))
    self.form = logform
    self.jac  = logjac
    self.pnms = pnames
    self.pmap = {}
    self.emap = {}
    self.parm = None
    self.perr = None
    self.cov  = None
  def _set_default_parms(self,xvals,yvals,evals):
    ''' Start by assuming the shift is zero. '''
    return (yvals[abs(xvals).argmin()],slope(np.log(xvals),yvals),0)

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
