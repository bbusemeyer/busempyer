import  mython as my
import subprocess as sp
import numpy as np
import os
import json

# Reads a qwalk input section.
def read_section(inp,key,pos):
  res = my.Ldict()
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
      print('Warning, reading commented lines is incomplete!')
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

def convert(s):
  return float(s)
  a=eval(s)
  m=float(a[0])
  for i in a:
    m=max(m,i)
  return m

def read_dm(inpf):
  """Read in the 1-RDM and/or 1-RDM and diagonals of the 2-RDM, if only 
  the diagonals were calculated"""
  nmo=0
  while True:
    line=inpf.readline()
    if line.find("tbdm")!=-1:
      spl=line.split()
      nmo=int(spl[2])
      break
    if line=="":
      return None
  #print "nmo ",nmo
  data={}
  #one-body up, one-body up err,   two-body-uu, two-body-uu err
  for nm in ['ou','oue','od','ode','tuu','tuue','tud','tude','tdu','tdue','tdd','tdde']:
    data[nm]=np.zeros((nmo,nmo))
  
  while True:
    line=inpf.readline()
    if line=="":
      break;
    if line.find("tbdm: states") != -1:
      data['states']=np.array(map(int,line.split()[3:-2]))
      #print data['states']
    if line.find("One-body density") != -1:
      line=inpf.readline()
      for i in range(0,nmo):
        for j in range(0,nmo):
          a=inpf.readline().split()
          data['ou'][i,j]=convert(a[2])
          data['oue'][i,j]=convert(a[3])
          data['od'][i,j]=convert(a[4])
          data['ode'][i,j]=convert(a[5])
    if line.find("two-body density") != -1:
      line=inpf.readline()
      #print "two-body",line
      for i in range(0,nmo):
        for j in range(0,nmo):
          a=inpf.readline().split()
          #print i,j,"a ",a
          data['tuu'][i,j]=convert(a[4])
          data['tuue'][i,j]=convert(a[5])
          data['tud'][i,j]=convert(a[6])
          data['tude'][i,j]=convert(a[7])
          data['tdu'][i,j]=convert(a[8])
          data['tdue'][i,j]=convert(a[9])
          data['tdd'][i,j]=convert(a[10])
          data['tdde'][i,j]=convert(a[11])
      break
  return data

def read_number_dens(inpf):
  data=np.zeros((0,0,0,0,0,0))
  data_err=np.zeros((0,0,0,0,0,0))
  while True:
    line=inpf.readline()
    #print line
    if line.find("Region fluctuation")!=-1:
      line=inpf.readline()
      spl=line.split()
      nspin=int(spl[1])
      maxn=int(spl[3])
      nregion=int(spl[5])
      data=np.zeros((nspin,nspin,nregion,nregion,maxn,maxn))
      #print "data ", data.shape
      data_err=np.zeros((nspin,nspin,nregion,nregion,maxn,maxn))
      
      for s1 in range(0,nspin):
        for s2 in range(0,nspin):
          for r1 in range(0,nregion):
            for r2 in range(0,nregion):
              line=inpf.readline()
              #print line
              for n1 in range(0,maxn):
                spl=inpf.readline().split()
                for n2 in range(0,maxn):
                  data[s1,s2,r1,r2,n1,n2]=float(spl[n2])
                for n2 in range(maxn,2*maxn):
                  data_err[s1,s2,r1,r2,n1,n2-maxn]=float(spl[n2])
      break
    elif line == '':
      return None,None
  return data,data_err

def moments(data,data_err):
  nspin=data.shape[0]
  nregion=data.shape[2]
  maxn=data.shape[4]
  
  avg=np.zeros((nspin,nregion))
  avg_err=np.zeros((nspin,nregion))
  for s in range(0,nspin):
    for r in range(0,nregion):
      for n in range(0,maxn):
        avg[s,r]+=n*data[s,s,r,r,n,n]
        avg_err[s,r]+=n**2*data_err[s,s,r,r,n,n]**2
  avg_err = avg_err**.5

  var=np.zeros((nspin,nregion))
  var_err=np.zeros((nspin,nregion))
  for s in range(0,nspin):
    for r in range(0,nregion):
      for n in range(0,maxn):
        var[s,r] += (n-avg[s,r])**2 * data[s,s,r,r,n,n]
        var_err[s,r]+=data_err[s,s,r,r,n,n]**2*(n-avg[s,r])**4 + 2*data[s,s,r,r,n,n]*avg_err[s,r]**2*(n-avg[s,r])**2
  var_err = var_err**.5

  covar=np.zeros((nspin,nspin,nregion,nregion))
  covar_err=np.zeros((nspin,nspin,nregion,nregion))
  for s1 in range(0,nspin):
    for s2 in range(0,nspin):
      for r1 in range(0,nregion):
        for r2 in range(0,nregion):
          corr=0.0
          corr_err=0.0
          for n1 in range(0,maxn):
            for n2 in range(0,maxn):
              p=data[s1,s2,r1,r2,n1,n2]
              pe=data_err[s1,s2,r1,r2,n1,n2]
              corr += p*(n2-avg[s2,r2])*(n1-avg[s1,r1])
              corr_err+= pe**2*(n1-avg[s1,r1])**2*(n2-avg[s2,r2])**2 + p**2*avg_err[s,r]**2*(n2-avg[s2,r2])**2 + p**2*avg_err[s,r]**2*(n1-avg[s1,r1])**2
          covar[s1,s2,r1,r2]=corr
          covar_err[s1,s2,r1,r2]=corr_err
  covar_err = covar_err**.5

  return avg,var,covar,avg_err,var_err,covar_err

# Use golsing to read out average energy and error.
def read_qenergy(logfile,gosling='./gosling'):
  statfilen = logfile.name.replace('.log','.stat')
  with open(statfilen,'w') as out:
    try:
      sp.call([gosling, logfile.name], stdout = out)
    except OSError:
      print("Cannot find gosling!")
      exit()
  with open(statfilen,'r') as F:
    for line in F:
      if 'total_energy0' in line:
        spl = line.split()
        return {'egy':spl[1],'err':spl[3],'var':spl[5]}
  print('ERROR: cannot find total_energy0 in stat file.')
  return {'egy':None,'err':None,'var':None}

def gen_optimize(dftfn,time='02:00:00'):
  loc = '/'.join([os.getcwd()]+dftfn.split('/')[:-1])
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

  return my.gen_qsub(exe,stdout=out,loc=loc,
                  name=loc+' gen_optimize',
                  time=time,
                  nn=2,np=12,
                  queue='secondary',
                  prep_commands=pc,
                  final_commands=fc)

#TODO add runtime checks for file dependence.
def gen_dmc(root,ts=0.01,jast=''):
  """ Generate DMC input file for a specified kpoint root. """

  # If no jast provided, will try to use one named like the dmc run.
  if jast == '':
    jast = root + '.opt.jast2'

  dmcfn = root+'.dmc'
  with open(dmcfn,'w') as dmcf:
    dmclines = []
    dmclines.append('method{ dmc ')
    dmclines.append('  timestep %s'%ts)
    dmclines.append('  save_trace %s.dmc.trace'%root)
    dmclines.append('}')
    dmclines.append(' include %s.sys'%root)
    dmclines.append(' trialfunc{ slater-jastrow ')
    dmclines.append(' wf1{ include %s.slater }'%root)
    dmclines.append(' wf2{ include %s }'%jast)
    dmclines.append('}')
    dmcf.write('\n'.join(dmclines))
  return dmcfn

#TODO make general to all systems.
# Copied from Lucas.
def gen_basis_orb(sysf,minorbf):
  basis_numbers={"FE":range(1,11),"SE":range(1,5)}

  totmonum=1
  atomnum=1
  for line in sysf:
    #assume that the ATOM section is all on one line
    #as is printed out by the converters
    if "ATOM" in line:
      spl=line.split()
      nm=spl[2]
      for b in basis_numbers[nm]:
        minorbf.write("%i %i %i 1\n"%(totmonum,b,atomnum))
        totmonum+=1
      atomnum+=1
  minorbf.write("COEFFICIENTS\n1.0\n\n")
  return totmonum

# Postprocess (compute density, fluctuations, and 1-RDM).
def gen_ppr(root,
            gosling='./gosling',
            jast='',
            minbasisfn='/home/busemey2/tools/mython/minbasis.json'):
  """ Generate postprocess input file for a specified kpoint root. """

  # If no jast provided, will try to use one named like the dmc run.
  if jast == '':
    jast = root + '.opt.jast2'

  pprfn = root+'.ppr'
  # To generate stat file.
  egy = read_qenergy(open(root+'.dmc.log','r'),gosling=gosling)
  with open(root+'.dmc.stat','r') as statf:
    for line in statf:
      if 'Threw' in line:
        nwarm = int(line.split()[4])

  # Output postprocess file.

  with open(pprfn,'w') as pprf:
    pprlines = []
    pprlines.append('method{ postprocess ')
    pprlines.append('readconfig %s.dmc.tracele'%root)
    pprlines.append('nskip %d'%(nwarm*2048))
    pprlines.append('density { density up   outputfile %s.dmc.up.cube }'%root)
    pprlines.append('density { density down outputfile %s.dmc.dn.cube }'%root)
    pprlines.append('average { region_fluctuation }')
    pprlines.append('average { tbdm_basis')
    pprlines.append('mode obdm')
    pprlines.append('include %s'%(root+'.min.basis'))
    pprlines.append('}')
    pprlines.append('}')
    pprlines.append('include %s.sys'%root)
    pprlines.append('trialfunc{ slater-jastrow ')
    pprlines.append('wf1{ include %s.slater }'%root)
    pprlines.append('wf2{ include %s }'%jast)
    pprlines.append('}')
    pprf.write('\n'.join(pprlines))

  # Generate minimum basis.
  with open(root+'.sys','r') as sysf, open(root+'.min.orb','w') as minorbf:
    nmo = gen_basis_orb(sysf,minorbf)
  with open(root+'.min.basis','w') as minbasisf:
    minbasis = json.load(open(minbasisfn,'r'))
    lines = []
    lines.append('orbitals {')
    lines.append('cutoff_mo')
    lines.append('orbfile %s.min.orb'%root)
    lines.append('nmo %d'%nmo)
    lines += minbasis['fese']
    lines.append('}')
    minbasisf.write('\n'.join(lines))

 
  return pprfn
