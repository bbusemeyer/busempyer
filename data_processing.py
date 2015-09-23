#!/usr/bin/python
<<<<<<< HEAD
import numpy as np
=======
from numpy            import array,linalg
from dm_tools         import read_dm
from read_numberfluct import read_number_dens,moments
from cubetools        import read_cube
from cryfiles_io      import read_cryinp, read_cryout, read_crytrace
from qfiles_io        import read_qfile, read_qenergy
from os               import getcwd
>>>>>>> 347f8307e8a313a9af804ef65d189a2bf78e8340
import os
import json
import cryfiles_io as cio
import qfiles_io as qio

# Temporary function to convert file names to metadata about the calculations.
# In the future, an optional metadata file should contain overall qualitiative
# descriptions of the data, like magnetic ordering, or pressure. These things
# are redundant, but convenient for interpreting sets of numbers together.
def convert_to_metadata(froot):
  mconv = {'che':'checkerboard',
           'str':'collinear',
           'unp':'unpolarized',
           'fer':'ferromagnetic',
           'bic':'bicollinear',
           'fst':'collinear, flip 1',
           'dim':'collinar, flip 2',
           'sta':'staggered',
           'cpx':'checkerboard',
           'spx':'collinear',
           'bpx':'bicollinear'}
  basename = froot.split('/')[-1]
  try:
    mag = mconv[basename[:3]]
  except:
    mag = None
  try:
    if 'px' in basename[:3]: # e.g. 'cpx'
      prs = float(basename.split('_')[1])
    else:
      prs = float(basename[3:].replace('_prs',''))
  except:
    prs = None
  return {'magnetic ordering':mag,'pressure':prs}

# Read a directory which has multiple files with data into a single dictionary
# with relevant information.
def read_dir(froot,gosling='./gosling'):
  """ Reads a CRYSTAL + QWalk directory's data into a dictionary.
  
  Current dictionary keys:
    fmixing, kdens, excited_energy_err, ordering, total_energy_err, 
    excited_energy, ts, tole, total_energy, se_height, tolinteg, 
    supercell, mixing, kpoint, spinlock, dft_energy, a, c, broyden, 
    dft_moments, 1rdm, access_root, average, covariance
  """

  ############################################################################
  # This first section should be edited to reflect naming conventions!
  dftfile   = froot+'.d12'
  dftoutf   = froot+'.d12.out'
  # Start by taking only the real k-points, since I'm sure these are on solid
  # ground, and have enough sample points. TODO generalize
  realk1 = array([1,3,8,10,27,29,34,36]) - 1
  realk2 = array([1,4,17,20,93,96,109,112]) - 1
  oldrealk = array(['k0','k1','k2','k3','k4','k5','k6','k7'])
  ############################################################################

  bres = {} # Data that is common to all k-points.
  ress = {} # Dict of all k-point data in directory.

  bres['access_root'] = os.getcwd() + '/' + froot

  print "Working on",froot+"..." 
  try:
    with open(froot+'_metadata.json','r') as metaf:
      metad = json.load(metaf)
    try:
      if metad['magnetic ordering'] != None:
        bres['ordering'] = metad['magnetic ordering']
    except KeyError: pass
    try:
      if metad['pressure'] != None:
        bres['pressure'] = metad['pressure']
    except KeyError: pass
    try:
      if metad['pressure'] != None:
        bres['pressure'] = metad['pressure']
    except KeyError: pass
  except IOError:
    print "  Didn't find any metadata"
  
  print "  DFT params and results..." 
  try:
    dftdat = cio.read_cryinp(open(dftfile,'r'))
    bres['a'] = dftdat['latparms'][0]
    bres['c'] = dftdat['latparms'][1]
    bres['se_height'] = dftdat['apos'][dftdat['atypes'].index(234)][-1]
    for key in ['mixing','broyden','fmixing','tolinteg',
                'kdens','spinlock','supercell','tole']:
      bres[key] = dftdat[key]
    dftdat = cio.read_cryout(open(dftoutf,'r'))
    bres['dft_energy'] = dftdat['dft_energy']
    bres['dft_moments'] = dftdat['dft_moments']
  except IOError:
    print "There's no dft in this directory!"
    return {}

  # Determine k-point set and naming convention.
  if os.path.isfile(froot+'_'+str(oldrealk[0])+'.sys'):
    realk=oldrealk
    print "Using old-style kpoint notation."
  elif os.path.isfile(froot+'_'+str(realk2[-1])+'.sys'):
    realk=realk2
    print "Using 6x6x6 kpoint notation."
  else:
    realk=realk1
    print "Using 4x4x4 kpoint notation."

  for rk in realk:
    kroot = froot + '_' + str(rk)

    print "  now DMC:",kroot+"..." 
    try:
      sysdat = qio.read_qfile(open(kroot+'.sys','r'))
      dmcinp = qio.read_qfile(open(kroot+'.dmc','r'))
      ress[rk] = bres.copy()
      ress[rk]['kpoint'] = sysdat['system']['kpoint']
      ress[rk]['ts'] = dmcinp['method']['timestep']
    except IOError:
      print "  (cannot find QMC input, skipping)"
      continue
    
    print "  energies..." 
    try:
      inpf = open(kroot+'.dmc.log','r')
      egydat = qio.read_qenergy(inpf,gosling)
      ress[rk]['dmc_energy']     =  egydat['egy']
      ress[rk]['dmc_energy_err'] =  egydat['err']
    except IOError:
      print "  (cannot find ground state energy log file)"

    try:
      inpf = open(kroot+'.ogp.log','r')
      ogpdat = qio.read_qenergy(inpf,gosling)
      ress[rk]['dmc_excited_energy']     =  ogpdat['egy']
      ress[rk]['dmc_excited_energy_err'] =  ogpdat['err']
    except IOError:
      print "  (cannot find excited state energy log file)"

    print "  densities..." 
    try:
      inpf = open(kroot+'.dmc.up.cube','r')
      upcube = read_number_dens(inpf)
      inpf = open(kroot+'.dmc.dn.cube','r')
      dncube = read_number_dens(inpf)
      ress[rk]['updens'] = upcube
      ress[rk]['dndens'] = dncube
    except IOError:
      print "  (cannot find electron density)"

    print "  fluctuations..." 
    try:
      inpf = open(kroot+'.ppr.o','r')
<<<<<<< HEAD
      fludat, fluerr = qio.read_number_dens(inpf)
      if fludat==None:
=======
      fludat, fluerr = read_number_dens(inpf)
      if fludat is None:
>>>>>>> 347f8307e8a313a9af804ef65d189a2bf78e8340
        print "  (Error in number fluctuation output, skipping)"
      else:
        avg, var, cov, avge, vare, cove = qio.moments(fludat,fluerr)
        ress[rk]['average']    = avg
        ress[rk]['covariance'] = cov
    except IOError:
      print "  (cannot find number fluctuation)"

    print "  1-RDM..." 
    try:
      inpf = open(kroot+'.ordm.o')
<<<<<<< HEAD
      odmdat = qio.read_dm(inpf)
      if odmdat==None:
=======
      odmdat = read_dm(inpf)
      if odmdat is None:
>>>>>>> 347f8307e8a313a9af804ef65d189a2bf78e8340
        print "  (Error in 1-RDM output, skipping)"
      else:
        ress[rk]['1rdm'] = odmdat
    except IOError:
      print "  (cannot find 1-RDM)"

    print "  done."; 

  return ress

def trace_analysis(dftfns,ids=[]):
  """
  Generate a dictionary of energy traces from DFT output, useful for checking
  how certain parameters affect convergence.
  """
  res = {}
  if ids == []: ids = dftfns
  for di,dftfn in enumerate(dftfns):
    with open(dftfn,'r') as dftf:
      res[ids[di]] = cio.read_crytrace(dftf)
  return res
