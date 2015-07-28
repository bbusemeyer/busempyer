import numpy as np

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
