# Read in CYRSTAL output file, output nearest neighbor list
import sys

finp = sys.argv[1]

while realine != '':
  line = finp.readline()
  if 'NEIGHBORS' not in line: continue
  for i in range(3): finp.readline()
  line = finp.readline().split()

  # CRYSTAL format is too annoying to continue!
