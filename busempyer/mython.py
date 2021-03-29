import subprocess as sub
import numpy as np
import scipy.linalg as lin
from json import JSONEncoder
from os import getcwd,rename
from datetime import datetime

# Lines should be a list of lists of words.
# Words are separated by spaces, lines by \n.
def lines2str(lines):
  outlines = []
  for line in lines:
    outlines.append(' '.join(map(str,line)))
  outstr = '\n'.join(outlines)
  return outstr

def cross_prod(sets):
  mg = np.array(np.meshgrid(*sets))
  mg = mg.swapaxes(0,-1)
  mg = mg.reshape(np.prod(mg.shape[:-1]),mg.shape[-1])
  return mg

def fix_duped_json(jsonfn,key='tbdm_basis',newfn=None):
  ''' When QWalk dumps a JSON with multiple keys from the same method, it just
  prints them several times. When a JSON reader reads that, it overwrites the
  previous keys. 

  As a hack fix, just go through and enumerate the repeated keys.'''
  if newfn is None: newfn='oiremcoidsajfdsa.json'

  newf=open(newfn,'w')
  keycount=0
  for line in open(jsonfn,'r'):
    if key in line:
      newf.write(line.replace(key,"%s_%d"%(key,keycount)))
      keycount+=1
    else:
      newf.write(line)

  newf.close()

  if newfn=='oiremcoidsajfdsa.json':
    rename(newfn,jsonfn)

class Ldict(dict):
  """
  Dictionary that can append data rather than overwriting it when append item is
  used more than once.
  """
  def __init__(self,**kwargs):
    dict.__init__(self,kwargs)

  def append(self,key,value):
    if key in dict.keys(self):
      if isinstance(dict.__getitem__(self,key),list):
        dict.__getitem__(self,key).append(value)
      elif isinstance(dict.__getitem__(self,key),bool):
        dict.__setitem__(self,key,value)
      else:
        dict.__setitem__(self,key,[dict.__getitem__(self,key),value])
    else:
      dict.__setitem__(self,key,value)

def format_nested_dict(nested,indent=""):
  """ Display all keys in a nested dictionary. 
  Don't use indent if calling on root."""
  outstr = ""
  for key in nested.keys():
    outstr += "{0}{1}\n".format(indent,key)
    if type(nested[key])==dict:
      outstr += format_nested_dict(nested[key],indent="  "+indent)
  return outstr

# Source (simplified by me):
# http://stackoverflow.com/questions/3488934/simplejson-and-numpy-array/24375113#24375113
class NumpyEncoder(JSONEncoder):
  def default(self, obj):
    if isinstance(obj, np.ndarray):
      print(obj)
      return dict(__ndarray__=obj.tolist(),
                  dtype=str(obj.dtype),
                  shape=obj.shape)
    # Let the base class default method raise the TypeError
    return JSONEncoder(self, obj)

class NumpyToListEncoder(JSONEncoder):
  def default(self, obj):
    if isinstance(obj, np.ndarray):
      return obj.tolist()
    # Let the base class default method raise the TypeError
    return JSONEncoder(self, obj)

def json_numpy_obj_hook(dct):
  """ Hook for facier numpy encoder

  Example usage:
  expected = np.arange(100, dtype=np.float)
  dumped = json.dumps(expected, cls=NumpyEncoder)
  result = json.loads(dumped, object_hook=json_numpy_obj_hook)
  """

  if isinstance(dct, dict) and '__ndarray__' in dct:
    data = np.ndarray(dct['__ndarray__'],dtype=dct['dtype'])
    return data.reshape(dct['shape'])
  return dct

def gen_qsub(exe,stdout='',loc='',name='',time='72:00:00',nn=1,np='allprocs',
    queue='batch', prep_commands=[],final_commands=[]):
  """ Generate a qsub file.
  
  Blank strings will generate useful defaults."""

  if stdout=='': stdout='stdout'
  if loc=='': loc=getcwd()
  if name=='': name=str(datetime.now()).replace(' ','_')
  header = []
  header.append('#!/bin/bash')
  if np=='allprocs': header.append('#PBS -l nodes=%d,flags=allprocs'%nn)
  else:              header.append('#PBS -l nodes=%d:ppn=%s'%(nn,str(np)))
  header.append('#PBS -q %s'%queue)
  header.append('#PBS -l walltime=%s'%time)
  header.append('#PBS -j oe')
  header.append('#PBS -m n')
  header.append('#PBS -N %s'%name)
  header.append('#PBS -o %s.qsub.out'%stdout)
  if np=='allprocs':
    exeline = 'mpirun %s &> %s'%(exe, stdout)
  elif nn*np > 1:
    exeline = 'mpirun -n %d %s &> %s'%(nn*np, exe, stdout)
  else:
    exeline = '%s &> %s'%(exe, stdout)
  commands = header + ['cd %s'%loc] + prep_commands + [exeline] + final_commands
  outstr = '\n'.join(commands)
  with open(loc+'/qsub.in','w') as qsin:
    qsin.write(outstr)
  return loc+'/qsub.in'

def gen_qsub_general(exelines,loc='',name='',time='72:00:00',nn=1,np='allprocs',
    queue='batch', prep_commands=[],final_commands=[]):
  """ Generate a qsub file.
  
  Blank strings will generate useful defaults. 'general' is a new version that
  allows different execution strings."""

  if loc=='': loc=getcwd()
  if name=='': name=str(datetime.now()).replace(' ','_')
  header = []
  header.append('#!/bin/bash')
  if np=='allprocs': header.append('#PBS -l nodes=%d,flags=allprocs'%nn)
  else:              header.append('#PBS -l nodes=%d:ppn=%s'%(nn,str(np)))
  header.append('#PBS -q %s'%queue)
  header.append('#PBS -l walltime=%s'%time)
  header.append('#PBS -j oe')
  header.append('#PBS -m n')
  header.append('#PBS -N %s'%name)
  header.append('#PBS -o %s/%s.qsub.out'%(loc,name))
  commands = header + ['cd %s'%loc] + prep_commands + exelines + final_commands
  outstr = '\n'.join(commands)
  with open(loc+'/qsub.in','w') as qsin:
    qsin.write(outstr)
  return loc+'/qsub.in'

class Bootstrapper:
  """ Can be used to resample some function over and over to gain information
  about the statistics. Currenly only computes the variance. """
  def __init__(self,func,resample):
    """ Definition of inputs:
    func: function that accepts random variable inputs.
    resampler: Generates samples of input from whatever distribution you deem fit. """
    self.func = func
    self.resample = resample
  def gen_stats(self,nsamples=1000):
    """ Compute statistics of func results from distribution of resampler. """
    stats = {}
    mom1 = np.array(self.func(self.resample()))
    mom2 = np.array(self.func(self.resample()))**2
    for i in range(nsamples-1):
      mom1 += np.array(self.func(self.resample()))
      mom2 += np.array(self.func(self.resample()))**2
    stats['mean'] = mom1/nsamples
    stats['variance'] = mom2/nsamples - stats['mean']**2
    return stats

  def gen_stats_parallel(self,nsamples=1000,ncores=8):
    """ Compute statistics of func results from distribution of resampler. 
    Computes in parallel, which is *sometimes* faster (for harder func calls). """
    raise AssertionError("Not implemented!")

def gaussian_matrix_resample(values,stdevs):
  return np.random.randn(*values.shape)*stdevs + values

class Bootstrapper_eigh(Bootstrapper):
  """ Bootstrapper for finding eigenvalues of hermitian matrices."""
  def __init__(self,matrix,error,overlap=None,overlap_err=None):
    if overlap_err is None:
      self.func = lambda mat:lin.eigh(mat,overlap)[0]
      self.resample = lambda: gaussian_matrix_resample(matrix,error)
    else:
      self.func = lambda mats:lin.eigh(mats[0],mats[1])[0]
      self.resample = lambda: (gaussian_matrix_resample(matrix,error),
                                 gaussian_matrix_resample(overlap,overlap_err))

class Bootstrapper_eig(Bootstrapper):
  """ Bootstrapper for finding eigenvalues of hermitian matrices."""
  def __init__(self,matrix,error,overlap=None,overlap_err=None):
    if overlap_err is None:
      self.func = lambda mat:lin.eig(mat,overlap)[0]
      self.resample = lambda: gaussian_matrix_resample(matrix,error)
    else:
      self.func = lambda mats:lin.eig(mats[0],mats[1])[0]
      self.resample = lambda: (gaussian_matrix_resample(matrix,error),
                                 gaussian_matrix_resample(overlap,overlap_err))

def run_output(command,outf=None):
  # Python 3.5
  #completed = sub.run(command,shell=True,
  #    stdout=sub.PIPE,stderr=sub.PIPE,universal_newlines=True)
  #if outf is not None:
  #  outf.write(completed.stdout)
  #  outf.write(completed.stderr)
  # Python 3.4
  completed = sub.check_output(command,shell=True,
      stderr=sub.STDOUT,universal_newlines=True)
  if outf is not None:
    outf.write(completed)
  return completed
