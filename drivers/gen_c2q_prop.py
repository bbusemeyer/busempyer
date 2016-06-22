import sys
import cryfiles_io as cryio

if len(sys.argv) < 2:
  print("Usage: gen_c2q_prop.py crystal_input [prop_outfile]")
  raise AssertionError

cryinp = cryio.read_cryinp(sys.argv[1])
kdens = cryinp['kdens']

outlines = [
    "NEWK",
    "{} {}".format(kdens,2*kdens),
    "CRYAPI_OUT",
    "END"
  ]

if len(sys.argv) == 3:
  with open(sys.argv[2],'w') as outf:
    outf.write('\n'.join(outlines))
else:
  print('\n'.join(outlines))
