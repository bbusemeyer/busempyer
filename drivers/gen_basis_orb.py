#!/usr/bin/python
# Original author: Lucas K. Wagner.
# Create an orb file from a sys file and a number of basis functions. Each
# is simply the corresponding basis element.
import sys

print("You'll need to input the element names in order, followed by the number present in the material.")
done = False
basis_numbers = {}
while not done:
  ion_name = input("Element name (empty for done):")
  if ion_name != "":
    ion_number = int(input("Number of %s ions:"%ion_name))
    basis_numbers[ion_name] = range(1,ion_number+1)
  else:
    done = True

# If you're reading this because you're annoyed at using the interface, go
# ahead and hard-code the dict here:
#basis_numbers={"FE":range(1,11),"SE":range(1,5)}

f=open(sys.argv[1],'r')

totmonum=1
atomnum=1
for line in f:
  #assume that the ATOM section is all on one line
  #as is printed out by the converters
  if "ATOM" in line:
    spl=line.split()
    nm=spl[2]
    for b in basis_numbers[nm]:
      print "%i %i %i 1"%(totmonum,b,atomnum)
      totmonum+=1
    atomnum+=1
print "COEFFICIENTS\n1.0\n\n"

