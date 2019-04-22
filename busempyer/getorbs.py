# Brian Busemeyer
# Extract up and down orbitals, write plot method for qwalk
from pickle import load

def main():
  conv = load(open('conv.pkl','rb'))
  out = from_conv(conv.orbitals[0],conv.system,'crys_0.orb',range(1,45))
  print(out)

def from_obj(orbitals,system,orbfile,orb_indices):
  ''' Make a plot input file using a ConverterManager. 
  Args: 
    orbital (Orbitals): see qwalk_objects/orbitals.py.
    system (System): see qwalk_objects/system.py.
    convmgr (str): path to converter pickle.
    orbfile (str): path to orbfile.
    orb_indices (array-like): orbitals to plot. NOTE THIS IS 1-BASED INDEXING!!
  Returns:
    str: input file for QWalk.
  '''
  orblines = orbitals.export_qwalk_orbitals(orbfile).split('\n')
  syslines = system.export_qwalk_sys().split('\n')

  outlines = [
      'method { plot',
      '  resolution 0.1',
      '  plotorbitals { ' + ' '.join(map(str,orb_indices)) + ' }'
    ] + ['  ' + line for line in orblines] + [
      '}'
    ] + syslines

  return '\n'.join(outlines)

if __name__ == "__main__":
  main()
