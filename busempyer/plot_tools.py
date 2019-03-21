''' Convenience tools for plotting.'''
import numpy as np
import os
from copy import deepcopy as copy
import matplotlib.pyplot as plt

### Matplotlib tools.

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

Favorite colors to use from pc:
  't':    '#1b9e77', # Teal.
  'do':   '#d95f02', # Dark orange.
  'g' :   '#66a61e', # Green.
  'pink': '#e7298a', # Dark pink.
  'dy':   '#e6ab02', # Dark yellow/gold.
  'dgray':'#666666', # Dark gray.
  'tan' :' #b15928', # Cowhide tan.
  'lp':   '#7570b3', # light purple.


Color maps:
  Good ``uniform'' sequential: viridis,plasma.
  Good diverging: Spectral,seismic,bwr,BrBG
"""

# Nicer colors for plotting from colorbrewer2.org.
# Common colors are first letter only, prefix with l or d means light or dark.
pc = {
    # 12-Class paired.
    'lb':'#a6cee3', # Light blue.
    'b' :'#1f78b4', # Blue.
    'lg':'#b2df8a', # Light green.
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
    'g' :   '#66a61e', # Green.
    'lp':   '#7570b3', # light purple.
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
    "D",
    "p",
    "d",
    "h",
    "v",
    "^",
    "<",
    ">",
    "H",
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

def gen_rgb(hexstr):
  ''' Make a tuple of RGB values from a hex string. '''
  hstr=hexstr.lstrip('#')
  return tuple(int(hstr[i:i+2], 16) for i in (0, 2 ,4))

def gen_hex(rgb):
  ''' Take an RGB tuple and convert to a #XXXXXX string.'''
  digits=[('0'+hex(rgb[i]).replace('0x',''))[-2:] for i in range(3)]
  return '#{}{}{}'.format(*digits)

def gen_spectrum(start,end,numpoints):
  ''' Sequential colors can be chosen automatically.'''
  start=np.array(gen_rgb(start))
  end=np.array(gen_rgb(end))
  dom=np.linspace(0,1,numpoints)[:,None]
  spectrum=(1-dom)*start[None,:] + dom*end[None,:]
  spectrum=spectrum.round().astype(int)
  spectrum=[gen_hex(spectrum[i,:]) for i in range(numpoints)]
  return spectrum

def make_plotargs(**kwargs):
  ''' Make a nice set of defaults for plot aesthetics.
  Args: 
    kwargs: additional arguments to modify defaults.
  '''
  # Defaults.
  myplotdef={
      'mew':0.5,
      'mec':'k',
      'ms':5,
      'lw':1
    }
  for arg in kwargs:
    myplotdef[arg] = kwargs[arg]
  return myplotdef

def make_errargs(**kwargs):
  ''' Make a nice set of defaults for plot aesthetics.
  Args: 
    kwargs: additional arguments to modify defaults.
  '''
  # Defaults.
  myerrdef={
      'capthick':1,
      'capsize':2,
      'ecolor':'k',
      'fmt':'none'
    }
  for arg in kwargs:
    myerrdef[arg] = kwargs[arg]
  return myerrdef

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

def fix_lims(ax_inp,factor=0.04,do_x=True,do_y=True):
  """
  Create buffer around all points that is factor*(data range) wide.
  Fixes it in-place.

  Args:
    ax_inp (matplotlib.Axes or Axes array): Axes to fix.
    factor (float): fraction of empty space around points.
    do_x (bool): Fix the x axis.
    do_y (bool): Fix the y axis.
  """

  minx,miny = np.Inf,np.Inf
  maxx,maxy = -np.Inf,-np.Inf
  if type(ax_inp) != np.ndarray:
    ax_array = np.array((ax_inp))
  else:
    ax_array = ax_inp
  xdata=np.array(())
  ydata=np.array(())
  for ax in ax_array.flatten():
    for line in ax.get_lines():
      # axvlines store the data as lists and often should be ignored.
      if type(line.get_xdata()) is type([]):
        continue
      xdata=np.concatenate((xdata,line.get_xdata()))
      ydata=np.concatenate((ydata,line.get_ydata()))

  maxx=xdata[~np.isnan(xdata)].max()
  maxy=ydata[~np.isnan(ydata)].max()
  minx=xdata[~np.isnan(xdata)].min()
  miny=ydata[~np.isnan(ydata)].min()

  xs = factor*(maxx-minx)
  ys = factor*(maxy-miny)
  for ax in ax_array.flatten():
    if xs != 0 and do_x: ax.set_xlim(minx-xs,maxx+xs)
    if ys != 0 and do_y: ax.set_ylim(miny-ys,maxy+ys)

def slope(x,y): return (y[-1]-y[0])/(x[-1]-x[0])

def fix_xticks(ax,**kwargs):
  ''' Convenience function for thin_ticks.'''
  ax.set_xticks(thin_ticks(ax.get_xticks(),kwargs))

def fix_yticks(ax,**kwargs):
  ''' Convenience function for thin_ticks.'''
  ax.set_yticks(thin_ticks(ax.get_yticks(),kwargs))

def fix_ticks(ax,**kwargs):
  ''' Convenience function for thin_ticks.'''
  fix_xticks(ax,kwargs)
  fix_yticks(ax,kwargs)

def thin_ticks(ticks,div=2,start=0,shift=0,append=0):
  ''' Remove from list at evenly spaced intervals.
  Args:
    div (int): Total number of ticks divided by this.
    start (int): Start the ticks from this point. 
    shift (int): Shift the starting tick by this much.
    append (int): Add this more ticks.
  Returns:
    list: new ticks.
  '''
  newticks = [ticks[div*i-shift] for i in range(start,len(ticks)//div+append)]
  return newticks

def idxmap(arraylike):
  """ Map elements of an array to it's index (reverse mapping of array itself)."""
  return dict(zip(arraylike,range(len(arraylike))))

def safemap(di,key):
  if key in di.keys():
    return di[key]
  else:
    return key

# The mother of all plotting tools.
class CategoryPlot: 
  def __init__(self,df,
      row='categoryplotdummy',col='categoryplotdummy',
      color='categoryplotdummy',mark='categoryplotdummy',
      labmap={},cmap=None,mmap=None,sharex=False,sharey=False):
    '''
    Use a pandas DataFrame to make plots broken down by color, row, column,
    and marker. Somewhat similar to what ggplot can handle (more elegantly).

    For example: df.columns=[x,y,z],
    cp=CategoryPlot(df,color='z',mark='z')
    cp.plot('x','y')

    Now cp.fig will have a figure of x vs y with color and marker broken down by z.
    
    Call self.plot() to actually make a plot. 
   
      Args:
        row: rows will differ by this quantity (default to one row).
        col: columns will differ by this quantity (default to one column).
        color: colors will differ by this quantity (default to one color).
        mark: markers will differ by this quantity (default to one marker).
        labmap: labels of data values are mapped using labmap first. Not in map means leave as-is.
        cmap: data values are mapped to these colors (default to ps['dark8']).
        mmap: data values are mapped to these markers (default to pm).
        sharex: x-axes are set to same limits.
        sharey: y-axes are set to same limits.
    '''



    if 'categoryplotdummy' in df.columns:
      print("CategoryPlot: Warning, I'm not going to use the 'categoryplotdummy' column!")


    assert df.shape[0]>0 and df.shape[1]>0, "Empty dataframe!"
    self.fulldf=df
    self.fulldf['categoryplotdummy']='categoryplotdummy'
    self.row=row
    self.col=col
    self.color=color
    self.mark=mark
    self.labmap=labmap
    self.plotargs={}

    if cmap is None:
      self.unique_colors=self.fulldf[color].unique()
      nc=self.unique_colors.shape[0]
      self.cmap=dict(zip(self.unique_colors,( (1+nc//(len(ps['dark8'])+len(ps['cb12'])))*(ps['dark8']+ps['cb12']) )[:nc]))
    else: 
      self.cmap=cmap
    self.cmap['categoryplotdummy']='none'
    
    if mmap is None:
      self.unique_marks=self.fulldf[mark].unique()
      nm=self.unique_marks.shape[0]
      self.mmap=dict(zip(self.unique_marks,pm[:self.unique_marks.shape[0]]))
      self.mmap=dict(zip(self.unique_marks,((1+nm//len(pm))*pm)[:nm]))
    else: 
      self.mmap=mmap
    self.mmap['categoryplotdummy']='s'

    self.labmap=lambda x:safemap(labmap,x)

    self.fig,self.axes=plt.subplots(
        self.fulldf[row].unique().shape[0],
        self.fulldf[col].unique().shape[0],
        squeeze=False,sharex=sharex,sharey=sharey
      )
    self.fig.set_size_inches(3*self.axes.shape[1]+0.5,
                         2.5*self.axes.shape[0]+0.5)
    self.rowmap=idxmap(df[row].unique())
    self.colmap=idxmap(df[col].unique())

  def plot(self,xvar,yvar,evar=None,plotargs={},errargs={},
      labrow=None,labcol=None,labloc=(0.7,0.9),
      fill=True,line=False,
      xscale='linear',yscale='linear'):
    '''
    Plot some freakin' data. Jeez how complicated is this object?!

    Args:
      xvar (str): column name for x-axis.
      yvar (str): column name for y-axis.
      evar (str): optional column name for errorbars.
      plotargs (dict): additonal options for matplotlib plot.
      errargs (dict): additional options for matplotlib errorbar.
      labrow (None or str): Label rows automatically; specify location as 'title', 'axes', or 'figure'.
      labcol (None or str): Label columns automatically; specify location as 'title', 'axes', or 'figure'.
      labloc (tuple): Location of annotation for labrow/labcol in axes fraction.
      fill (bool): whether points are filled or empty of color.
      line (bool): whether to draw a line between all points.
      xscale (str): 'linear' or 'log'; scale of the x axis.
      yscale (str): 'lienar' or 'log'; scale of the y axis.
    '''

    self.plotargs=plotargs
    for lab,axdf in self.fulldf.groupby([self.row,self.col],sort=False):
      row,col=lab
      ax=self.axes[self.rowmap[row],self.colmap[col]]

      # This will handle work pertaining to a single Axis.
      self.subplot(ax,xvar,yvar,evar,axdf,plotargs,errargs,fill,line,xscale,yscale)

      # Handle locations of labels for row and col variables.
      labtitle = []
      labannotate = []
      if labrow=='axes':
        self.axes[self.rowmap[row],0].set_ylabel(self.labmap(row))
      elif labrow=='title': 
        labtitle.append("{}: {}".format(self.row,self.labmap(row)))
      elif labrow=='figure':
        labannotate.append("{}: {}".format(self.row,self.labmap(row)))
      if labcol=='axes':
        self.axes[-1,self.colmap[col]].set_xlabel(self.labmap(col))
      elif labcol=='title': 
        labtitle.append("{}: {}".format(self.col,self.labmap(col)))
      elif labcol=='figure':
        labannotate.append("{}: {}".format(self.col,self.labmap(col)))
      if len(labtitle):
        ax.set_title('\n'.join(labtitle))
      if len(labannotate):
        print('\n'.join(labannotate))
        ax.annotate('\n'.join(labannotate),labloc,xycoords='axes fraction')

      # I'm a 90's baby.
      self.fig.tight_layout()

  def subplot(self,ax,xvar,yvar,evar=None,axdf=None,plotargs={},errargs={},
      fill=True,line=False,xscale='linear',yscale='linear'):
    ''' See plot. args are the same, but for only one plot in the grid.
    Additonal Args:
    ax (Axes): Axes instance to make a plot on.
    '''

    if xscale=='linear':
      if yscale=='linear':
        method=ax.plot
      elif yscale=='log':
        method=ax.semilogy
      else: raise ValueError("yscale should be linear or log, not %s"%yscale)
    elif xscale=='log':
      if yscale=='linear':
        method=ax.semilogx
      elif yscale=='log':
        method=ax.loglog
      else: raise ValueError("yscale should be linear or log, not %s"%yscale)
    else: raise ValueError("xscale should be linear or log, not %s"%xscale)

    self.plotargs=plotargs
    if axdf is None: axdf=self.fulldf
    for lab,df in axdf.groupby([self.mark,self.color],sort=False):
      mark,color=lab

      # Handle missing marks and colors.
      if mark not in self.mmap:
        print("Warning: %s has no mark. Assigning it default '.'"%mark)
        self.mmap[mark]='.'
      if color not in self.cmap:
        print("Warning: %s has no color. Assigning it default 'k'"%color)
        self.cmap[color]='k'

      if line:
        method(df[xvar],df[yvar],'-',
            color=self.cmap[color],**plotargs)

      if evar is not None:
        ax.errorbar(df[xvar],df[yvar],df[evar],**errargs)

      if fill:
        method(df[xvar],df[yvar],self.mmap[mark],
            color=self.cmap[color],**plotargs)
      else:
        if 'mew' not in self.plotargs: self.plotargs['mew']=1
        if 'mec' in self.plotargs: save=self.plotargs.pop('mec')
        method(df[xvar],df[yvar],self.mmap[mark],
            color='none',
            mec=self.cmap[color],**self.plotargs)
        if 'mec' in self.plotargs: self.plotargs['mec']=save

  def add_legend(self,ax=None,labmap={},args={}):
    """ Make a legend for the markers and/or colors. labmap maps data to
    pretty labels. locargs is passed to axes.legend(). Returns prox for legend
    handles. If there are two legends, the args should be a list.
    Args:
      ax (Axes, tuple, or None): Different options:
        Axes--use this Axes instance to place the legend. 
        tuple--use self.axes[ax] to place the legend.
        None--use self.axes[0,0] to place the legend.
      labmap (dict): map categories to labels that will appear in legend.
      args (dict): other options for matplotlib's legend call. 
        If thre are two legends (when mark and color are different descrimiators for the data),
        this should be a tuple, one for mark and one for color).
    """
    if ax is None: 
      ax=self.axes[0,0]
    elif type(ax) == tuple:
      ax=self.axes[ax]


    # Legend's marks should be fully visible.
    safeargs=copy(self.plotargs)
    if 'mew' in self.plotargs:
      safeargs.pop('mew')

    legargs = copy(self.plotargs)
    legargs['alpha'] = 1.0

    # Minimize needed labels:
    if self.mark=='categoryplotdummy': 
      prox=[plt.Line2D([],[],
            linestyle='',
            marker=self.mmap['categoryplotdummy'],color=self.cmap[unique],label=self.labmap(unique),
            **legargs
          ) for unique in self.unique_colors
        ]
      leg=ax.legend(handles=prox,**args)
    elif self.color=='categoryplotdummy': 
      prox=[plt.Line2D([],[],
            linestyle='',
            marker=self.mmap[unique],color=self.cmap['categoryplotdummy'],label=self.labmap(unique),
            **legargs
          ) for unique in self.unique_marks
        ]
      leg=ax.legend(handles=prox,**args)
    elif self.color==self.mark:
      prox=[plt.Line2D([],[],
            linestyle='',
            marker=self.mmap[unique],color=self.cmap[unique],label=self.labmap(unique),
            **legargs
          ) for unique in self.unique_colors
        ]
      leg=ax.legend(handles=prox,**args)
    else:
      if type(args)==dict:
        args=[args,args]
      cprox=[plt.Line2D([],[],
            linestyle='',
            marker=self.mmap['categoryplotdummy'],color=self.cmap[unique],label=self.labmap(unique),
            **legargs
          ) for unique in self.unique_colors
        ]
      mprox=[plt.Line2D([],[],
            linestyle='',
            marker=self.mmap[unique],color=self.cmap['categoryplotdummy'],label=self.labmap(unique),mew=1.0,
            **safeargs
          ) for unique in self.unique_marks
        ]
      prox=cprox,mprox
      mlegend=ax.legend(handles=mprox,**(args[1]))
      mlegend._legend_box.align='left'
      ax.add_artist(mlegend)
      leg=[ax.legend(handles=cprox,**args[0]),mlegend]
      leg[0]._legend_box.align='left'
    return leg

  def axvline(self,*args,**kwargs):
    ''' Add a vline to all axes in the set.'''
    for ax in self.axes.ravel():
      ax.axvline(*args,**kwargs)

  def axhline(self,*args,**kwargs):
    ''' Add a hline to all axes in the set.'''
    for ax in self.axes.ravel():
      ax.axhline(*args,**kwargs)

  def set_xlabel(self,label,**kwargs):
    ''' Label all bottom axes.
    Args:
      label (str): what do you think?!
      kwargs: options to Axes.set_xlabel.
    '''
    for ax in self.axes[-1,:]:
      ax.set_xlabel(label,kwargs)

  def set_ylabel(self,label,**kwargs):
    ''' Label all left axes.
    Args:
      label (str): what do you think?!
      kwargs: options to Axes.set_xlabel.
    '''
    for ax in self.axes[:,0]:
      ax.set_ylabel(label,kwargs)

