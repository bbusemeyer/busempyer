import numpy as np
from scipy.optimize import curve_fit

def fix_lims(ax_array,factor=0.04):
  """
  Create buffer around all points that is factor*(data range) wide.
  """
  minx,miny = np.Inf,np.Inf
  maxx,maxy = -np.Inf,-np.Inf
  for ax in ax_array.flatten():
    for line in ax.get_lines():
      if line.get_data()[0].max() > maxx:
        maxx = line.get_data()[0].max()
      if line.get_data()[1].max() > maxy:
        maxy = line.get_data()[1].max()
      if line.get_data()[0].min() < minx:
        minx = line.get_data()[0].min()
      if line.get_data()[1].min() < miny:
        miny = line.get_data()[1].min()
  xs = factor*(maxx-minx)
  ys = factor*(maxy-miny)
  for ax in ax_array.flatten():
    if xs != 0: ax.set_xlim(minx-xs,maxx+xs)
    if ys != 0: ax.set_ylim(miny-ys,maxy+ys)

class FitFunc:
  """
  Define a function that can be fit to and plotted with.
  """
  def __init__(self,form,jacobian=None,pnames=[]):
    self.form = form
    self.jac  = jacobian
    self.pnms = pnames
    self.parm = None
    self.perr = None
    self.cov  = None

  def fit(self,df,xcol,ycol,ecol,*p0):
    """
    Use cols from df to fit params with initial values p0
    """
    fit = curve_fit(self.form,
      df[xcol],df[ycol],sigma=df[ecol],
      absolute_sigma=True,
      p0=p0)
    self.parm = np.array(p0)
    self.perr = np.array(p0)
    for pi,p in enumerate(p0):
      self.parm[pi] = fit[0][pi]
      self.perr[pi] = fit[1][pi][pi]**.5
    self.cov  = fit[1]

  def eval(self,x):
    """
    Evaluate fitted function at point x.
    """
    if self.parm == None: return None
    else:
      return self.form(x,*self.parm)
  def eval_error(self,x):
    """
    Error from evaluating fitted function at point x.
    """
    if self.parm==None or self.perr==None or self.jac==None: return None
    else:
      return np.dot( self.jac(x,*self.parm).T,
                     np.dot(self.cov,
                            self.jac(x,*self.parm)))**.5

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
    self.parm = None
    self.perr = None
    self.cov  = None

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
    self.parm = None
    self.perr = None
    self.cov  = None
