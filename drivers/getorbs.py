# Brian Busemeyer
# Extract up and down orbitals, write plot method for qwalk

def produce_inputs(fileroot):
  openf = open(fileroot+'.slater','r')
  upoutpf = open('orbplotup','w')
  dnoutpf = open('orbplotdn','w')
  uporbs = []
  dnorbs = []

  upoutpf.write('method { PLOT\n')
  dnoutpf.write('method { PLOT\n')
  upoutpf.write('  resolution 0.1 \n')
  dnoutpf.write('  resolution 0.1 \n')

  ln = openf.readline().split()[0]

  while ln.split()[0] != '}': 
    ln = openf.readline()
    upoutpf.write(ln)
    dnoutpf.write(ln)

  while ln != ['#Spin','up','orbitals']: 
    ln = openf.readline().split()
      
  ln = openf.readline().split()
  while ln != ['#Spin','down','orbitals']:
    uporbs += ln
    ln = openf.readline().split()

  ln = openf.readline().split()
  while ln != []:
    if ln[-1] == '}': del ln[-1]
    dnorbs += ln
    ln = openf.readline().split()

  upoutpf.write('PLOTORBITALS{ ')
  dnoutpf.write('PLOTORBITALS{ ')
  for oi in uporbs: upoutpf.write(oi+' ')
  for oi in dnorbs: dnoutpf.write(oi+' ')
  upoutpf.write('}')
  dnoutpf.write('}')
  upoutpf.write('\n}\n')
  dnoutpf.write('\n}\n')
  upoutpf.write('include '+fileroot+'.sys')
  dnoutpf.write('include '+fileroot+'.sys')

def produce_inputs_spinless(fileroot):
  openf = open(fileroot+'.slater','r')
  outpf = open('orbplot','w')
  orbs = []

  outpf.write('method { PLOT\n')
  outpf.write('  resolution 0.1 \n')

  ln = openf.readline().split()[0]

  while ln.split()[0] != '}': 
    ln = openf.readline()
    outpf.write(ln)

  while ln != ['#Spin','up','orbitals']: 
    ln = openf.readline().split()
      
  ln = openf.readline().split()
  while ln != ['#Spin','down','orbitals']:
    orbs += ln
    ln = openf.readline().split()

  ln = openf.readline().split()
  while ln != []:
    if ln[-1] == '}': del ln[-1]
    orbs += ln
    ln = openf.readline().split()

  outpf.write('PLOTORBITALS{ ')
  for oi in orbs: outpf.write(oi+' ')
  outpf.write('}')
  outpf.write('\n}\n')
  outpf.write('include '+fileroot+'.sys')

if __name__ == "__main__":
  from sys import argv
  if len(argv) < 2:
    raise AssertionError("""
    Not enough arguements.
    Usage: getorbs.py qwalk_root [spinless]""")
  fileroot = argv[1]
  if len(argv)>2 and argv[2]=="spinless":
    raise NotImplementedError # TODO
    produce_inputs_spinless(fileroot)
  else:
    produce_inputs(fileroot)
