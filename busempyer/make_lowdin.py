from pickle import load

def make_lowdin():
  conv = load(open('conv.pkl','rb'))
  orbs = conv.orbitals[0]
  sys = conv.system

  #orbs.basis['H'] = [b for b in orbs.basis['H'] if b['angular'].lower() in ('s','p')]
  orbs.basis['H'] = orbs.basis['H'][:2]

  orblines = make_minorb(orbs,sys)
  orbfn = 'lowdin_in.orb'
  nmo = len(orblines)-2
  with open(orbfn,'w') as outf:
    outf.write('\n'.join(orblines))

  outlines = [
      'method { lowdin',
      '  resolution 0.05',
      '  orbitals { ',
      '    blas_mo',
      '    orbfile {}'.format(orbfn),
      '    nmo {}'.format(nmo),
    ] + ['    ' + line for line in orbs.export_qwalk_basis().split('\n')] + [
      '  }',
      '}',
    ] + [sys.export_qwalk_sys()]

  return outlines

def make_minorb(orbs,sys):

  outlines = []
  bmap = {'s':1,'p':3,'5d':5}
  mocount = 1
  aocount = 1
  cecount = 1
  for atom in sys.positions:
    aocount = 1
    for btype in orbs.basis[atom['species']]:
      for lx in range(bmap[btype['angular'].lower()]):
        outlines += ['{} {} {} 1'.format(mocount,aocount,cecount)]
        mocount += 1
        aocount += 1
    cecount += 1
  outlines += ['COEFFICIENTS','1.0']

  return outlines

if __name__=='__main__':
  with open('lowdin','w') as outf:
    outf.write('\n'.join(make_lowdin()))
