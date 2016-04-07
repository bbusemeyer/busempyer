#!/usr/bin/python
# Original author: Lucas K. Wagner.
# Create an orb file from a sys file and a number of basis functions. Each
# is simply the corresponding basis element.
import sys

if len(sys.argv) != 4:
  print("Usage: ")
  print("python gen_basis_orb.py your.sys your.min.basis output.orb")

sysf=open(sys.argv[1],'r')
basfn=sys.argv[2]
outf=open(sys.argv[3],'w')

# Convert basis type to range of basis elements numbers.
basis_converter = {"S":1,"P":3,"5D":5}
def find_basis_numbers(element,basfn):
  basf = open(basfn,'r')
  num_basis = 0
  baswords = basf.read().split()
  elmflg = False
  basflg = False
  for wnum,word in enumerate(baswords):
    if word == element: elmflg = True
    if (word == "{") and (elmflg): basflg = True
    if elmflg and basflg and (word in basis_converter.keys()):
      num_basis+=basis_converter[word]
    if elmflg and basflg and word == "}":
      basf.close()
      return range(1,num_basis+1)
  basf.close()
  return "Error in find_num_basis!"


# Collect list of atoms from sys file.
atomlist = []
for line in sysf:
  #assume that the ATOM section is all on one line
  #as is printed out by the converters
  if "ATOM" in line:
    spl=line.split()
    atomlist.append(spl[2])

# Find number of basis elements for each atom.
basis_numbers = {}
for atom in set(atomlist):
  basis_numbers[atom] = find_basis_numbers(atom,basfn)

# Output orb file.
totmonum=1
atomnum=1
for atom in atomlist:
  for b in basis_numbers[atom]:
    outf.write("%i %i %i 1\n"%(totmonum,b,atomnum))
    totmonum+=1
  atomnum+=1
outf.write("COEFFICIENTS\n1.0\n")
