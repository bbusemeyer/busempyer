#!/usr/bin/env python3

from math import ceil
from numpy import asarray,arange,loadtxt,array_split
from busempyer.plot_tools import matplotlib_header,make_plotargs as pargs, pc
from matplotlib.pyplot import subplots
from h5py import File
from busempyer.afqmc_tools import read_raw_afqmc
from busempyer.qmcdata import estimate_warmup, reblock
from os.path import exists
from qharv.reel.forlib.stats import corr as compute_autocorr
from ase.units import Ha
import logging
matplotlib_header()

def interface():
  logging.basicConfig(level=logging.INFO)
  ''' Will try to figure out which code you're trying to plot based on the information given.'''
  from argparse import ArgumentParser
  ap = ArgumentParser('plotrace.py')
  ap.add_argument('--tracedata','-t',dest='tracedata',type=str,default=None,
      help="File name of data where energy trace is stored.")
  ap.add_argument('--name','-n',type=str,default="energytotal",
      help="Name of data in file, if important.")
  ap.add_argument('--drop_outliers','-d',dest='drop_outliers',action='store_true',
      help="Hide outliars in plot (just count and report the number of them instead).")
  ap.add_argument('--warmup','-w',dest='warmup',type=int,default=None,
      help="Manually input the warmup (otherwise will be estimated).")
  ap.add_argument('--preblock','-p',dest='preblock',type=int,default=0,
      help="Specify max initial number of blocks to reduce number of points to reblock.")
  ap.add_argument('--figname','-f',dest='figname',default='trace',
      help="Figure name to save pdf to.")
  args = ap.parse_args()

  # Case AFQMCLab
  if args.tracedata is None:
    assert exists('HNum.dat'), "Use -t to specify data, or else this code looks for HNum and den.dat for AFQMCLab."
    data,_ = read_raw_afqmc()
    trace = data['energy']
    itime = data['beta']
  # Case PyQMC.
  elif 'h5' in args.tracedata: 
    hdf = File(args.tracedata,'r')
    trace = hdf[args.name][()]
    itime = None
  else:
    trace = loadtxt(args.tracedata)
    itime = None

  plotrace(trace,itime,
      warmup=args.warmup,
      drop_outliers=args.drop_outliers,
      preblock=args.preblock,
      figname=args.figname,
    )

def plotrace(trace,itime=None,warmup=None,drop_outliers=False,preblock=1,figname="trace"):
  trace = asarray(trace)

  if 1 < preblock < trace.shape[0]:
    trace = asarray([a.mean() for a in array_split(trace, preblock)])

  if warmup is None:
    warmup = estimate_warmup(trace)
    if warmup < 0:
      print("** Warmup impossible to estimate. I'm setting warmup to 0 for debugging purposes. **")

  if itime is None:
    itime = arange(trace.shape[0])
    warmtime = warmup
  else:
    itime = asarray(itime)
    warmtime = warmup * (itime[1]-itime[0])

  ekeep  = trace[warmup:]
  itkeep = itime[warmup:]

  autotime = compute_autocorr(ekeep)
  blocks = reblock(ekeep, ceil(autotime))

  if drop_outliers:
    outliers = abs(ekeep - trace.mean()) > 10*trace.std()
    print(f"...there are {outliers.sum()} outliers!! They will be dropped from the plot...")
    itime = itime[~outliers]
    trace = trace[~outliers]

  #blockitime = int(round(itkeep.shape[0]/blocks.shape[0]))
  blockitime = itkeep[::ceil(autotime)][:blocks.shape[0]]

  mean = ekeep.mean()
  std  = ekeep.std(ddof=1)
  serr = std * (autotime / ekeep.shape[0])**0.5

  fig,ax = subplots(1,2,gridspec_kw={'width_ratios': [3, 1]},sharey=True)
  ax[0].plot(itime,trace,'.',color=pc['g'],alpha=0.3,**pargs())
  #print(blockitime, blocks)
  ax[0].plot(blockitime,blocks,'.',color=pc['r'],**pargs())
  ax[0].axvline(warmtime,color=pc['grey'],lw=1)
  ax[0].axhline(mean,color=pc['grey'],lw=1)
  #for mult in range(-3,4):
    #ax[0].axhline(mean + mult*serr,ls='--',color=pc['grey'],lw=0.5)
  
  ax[1].axhline(mean,color=pc['grey'],lw=1)
  ax[1].hist(ekeep,bins=50,orientation='horizontal',color=pc['g'])


  print(f"\nSTD: {std:0.6f} Ha, STERR: {serr:0.6f} Ha, Autocorrelation: {autotime}")
  print(f"     {std*Ha:0.5f}  eV         {serr*Ha:0.5f}  eV")

  print(f"Saving {figname}")
  fig.set_size_inches(4,3)
  fig.tight_layout()
  fig.savefig(f"{figname}.pdf")
 
if __name__=='__main__': interface()
