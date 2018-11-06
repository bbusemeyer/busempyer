#!/usr/bin/python
# Brian Busemeyer.
# Create table of structure data from a cif file for plotting and comparison.
from pymatgen.io.cifio import CifParser
import sys

if len(sys.argv) < 3:
  print "Usage:"
  print "{fn} <list-of-cifs>".format(fn=sys.argv[0])
  exit()

cifnames = sys.argv[1:]
outfn    = raw_input("Output datafile (default=struct.dat):")
if outfn == "": outfn = "struct.dat"
outlines = [' '.join(('cif','prop','val'))]
for cifname in cifnames:
  print cifname
  parser = CifParser(cifname)
  awkcif = parser.as_dict()
  if len(awkcif.keys()) > 1:
    print "More than one structure in this cif. Not implemented."
    exit()
  cifd   = awkcif[awkcif.keys()[0]]

  for latc in ['a','b','c']:
    key = '_cell_length_'+latc
    outlines += [' '.join((cifname,latc,cifd[key]))]

  for si,site in enumerate(cifd['_atom_site_label']):
    outlines+=[' '.join((cifname,
                         site.replace('1',''),
                         cifd['_atom_site_fract_z'][si]))]
  outlines+=[' '.join((cifname,'author',
    cifd['_publ_author_name'][0].split()[0].replace(',','')))]
  outlines+=[' '.join((cifname,'group',
    cifd['_symmetry_space_group_name_H-M'].replace(' ','')))]
  try: 
    outlines+=[' '.join((cifname,'pressure',
      cifd['_cell_measurement_pressure']))]
  except KeyError:
    outlines+=[' '.join((cifname,'pressure', '0.0'))]

with open(outfn,'w') as outf:
  outf.write('\n'.join(outlines))
