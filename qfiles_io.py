from mython           import Ldict
from subprocess       import call
from os               import getcwd
from mython           import gen_qsub

# Reads a qwalk input section.
def read_section(inp,key,pos):
  res = Ldict()
  while inp[pos] != '}':
    if isinstance(inp[pos],float):           # If it's a number.
      while isinstance(inp[pos],float):
        #print 'Found data',key,inp[pos]
        res.append(key,inp[pos])
        pos += 1
    elif inp[pos] == '{':                    # Else if it's demarking a section,
      if isinstance(inp[pos+1],str):         # which could have keywords,
        label = inp[pos+1].lower()
        pos += 2
        #print 'Reading section',key,label
        val,pos = read_section(inp,key,pos)
        if label != False:
          val['label'] = label
        res.append(key,val)
      else:                                  # or just numbers.
        pos += 1
        val = []
        while isinstance(inp[pos],float):
          val.append(inp[pos])
          pos += 1
        if len(val) == 1: val = val[0]
        #print 'Found data',key,val
        res[key] = val
        pos += 1
    else:                                    # Else it's a keyword.
      key = inp[pos].lower()
      if key not in res.keys():
        #print 'Setting',key
        res[key] = True
      pos += 1
  pos += 1
  return res,pos

# Reads a qwalk input file.
def read_qfile(inpf):
  inpstr = ''
  for line in inpf:
    if '#' in line: # TODO: Needs to be fixed when '#' isn't the first thing in line.
      print 'Warning, reading commented lines is incomplete!'
      continue
    inpstr += line
  # Ensure correct splitting. This is inefficient for large files.
  inpstr = inpstr.replace('{',' { ')
  inpstr = inpstr.replace('}',' } ')
  inp    = inpstr.split() + ['}'] # Now I can use "read_section" on the file!
  for i in range(len(inp)): # Make everything numeric into floats.
    try:                inp[i] = float(inp[i])
    except ValueError:  pass
  return read_section(inp,inpf.name,0)[0]

# Use golsing to read out average energy and error.
def read_qenergy(logfile,gosling='./gosling'):
  statfilen = logfile.name.replace('.log','.stat')
  with open(statfilen,'w') as out:
    try:
      call([gosling, logfile.name], stdout = out)
    except OSError:
      print "Cannot find gosling!"
      exit()
  with open(statfilen,'r') as F:
    for line in F:
      if 'total_energy0' in line:
        spl = line.split()
        return {'egy':spl[1],'err':spl[3],'var':spl[5]}
  print 'ERROR: cannot find total_energy0 in stat file.'
  return {'egy':None,'err':None,'var':None}

def gen_optimize(dftfn):
  loc = '/'.join([getcwd()]+dftfn.split('/')[:-1])
  root = dftfn.split('/')[-1].replace('.d12','')

  opt1fn = root+'_0.opt1'
  with open('/'.join((loc,opt1fn)),'w') as opt1f:
    opt1lines = []
    opt1lines.append('method{ optimize }')
    opt1lines.append('include %s_0.sys'%root)
    opt1lines.append('trialfunc{ slater-jastrow ')
    opt1lines.append('wf1{ include %s_0.slater }'%root)
    opt1lines.append('wf2{ include %s.jast2  }'%root)
    opt1lines.append('}')
    opt1f.write('\n'.join(opt1lines))

  optfn  = root+'_0.opt'
  with open('/'.join((loc,optfn)),'w') as optf:
    optlines = []
    optlines.append('method{ optimize }')
    optlines.append('include %s_0.sys'%root)
    optlines.append('trialfunc{ include %s_0.opt.wfout }'%root)
    optf.write('\n'.join(optlines))

  exe = '~/bin/qwalk $inpfile'
  out = optfn+'.out'
  pc =  ['module load openmpi/1.4-gcc+ifort']
  pc += ['if [ -f {root}_0.opt.wfout ]'.format(root=root)]
  pc += ['then inpfile={0}'.format(optfn)]
  pc += ['else inpfile={0}'.format(opt1fn)]
  pc += ['fi']
  fc = []
  fc += ['if [ ! -f {root}_0.opt.wfout ]'.format(root=root)]
  fc += ['then cp {root}_0.opt1.wfout {root}_0.opt.wfout'.format(root=root)]
  fc += ['fi']
  fc += ['~/bin/separate_jastrow {root}_0.opt.wfout > {root}_0.opt.jast2'.format(root=root)]

  return gen_qsub(exe,stdout=out,loc=loc,
                  name=loc+' gen_optimize',
                  time='02:00:00',
                  nn=2,np=12,
                  queue='secondary',
                  prep_commands=pc,
                  final_commands=fc)

#TODO add runtime checks for file dependence.
def gen_dmc(syslist,time=120,nproc=2048):
  locs = ['/'.join([getcwd()]+fn.split('/')[:-1]) for fn in syslist]
  roots = [fn.split('/')[-1].replace('.sys','') for fn in syslist]

  #TODO generalize when not gamma optimization.
  gamma = 'not found'
  for root in roots:
    if root.endswith('_0'):
      gamma = root
  if gamma == 'not found':
    print "Can't find gamma point calculation!"
    exit()

  info = zip(locs,roots)
  for (loc,root) in info:
    dmcfn = root+'.dmc'
    with open('/'.join((loc,dmcfn)),'w') as dmcf:
      dmclines = []
      dmclines.append('method{ dmc ')
      dmclines.append('timestep 0.01')
      dmclines.append('save_trace %s.dmc.trace'%root)
      dmclines.append('}')
      dmclines.append('include %s.sys'%root)
      dmclines.append('trialfunc{ slater-jastrow ')
      dmclines.append('wf1{ include %s.slater }'%root)
      dmclines.append('wf2{ include %s.opt.jast2 }'%gamma)
      dmclines.append('}')
      dmcf.write('\n'.join(dmclines))

  qsub = []
  qsub.append('qsub')
  qsub.append('-q prod')
  qsub.append('-A SuperMatSim')
  qsub.append('-t {time}'.format(time=time))
  qsub.append('-n {nproc}'.format(nproc=nproc))
  qsub.append('--mode c32')
  qsub.append('-o {gamma}.dmc.out'.format(gamma='/'.join((loc,gamma))))
  qsub.append('~/bin/qwalk')
  for loc,root in info:
    qsub.append(loc+'/'+root+'.dmc')
  return ' '.join(qsub)
