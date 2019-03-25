# Brian Busemeyer
# Extract up and down orbitals, write plot method for qwalk
from pickle import load

def main():
  out = from_conv('conv.pkl','crys_0.orb',range(1,45))
  print(out)

def from_conv(convmgr,orbfile,orbitals):
  ''' Make a plot input file using a ConverterManager. 
  Args: 
    convmgr (str): path to converter pickle.
    orbfile (str): path to orbfile.
    orbitals (array-like): orbitals to plot. NOTE THIS IS 1-BASED INDEXING!!
  Returns:
    str: input file for QWalk.
  '''
  conv = load(open(convmgr,'rb'))
  orblines = conv.orbitals[0].export_qwalk_orbitals(orbfile).split('\n')
  syslines = conv.system.export_qwalk_sys().split('\n')

  outlines = [
      'method { plot',
      '  resolution 0.1',
      '  plotorbitals { ' + ' '.join(map(str,orbitals)) + '}'
    ] + ['  ' + line for line in orblines] + [
      '}'
    ] + syslines

  return '\n'.join(outlines)

if __name__ == "__main__":
  main()
