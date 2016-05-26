import sys
import os

# Takes a list of autogen directories, and prints out which of them have
# corrupted DMC trace files.

if len(sys.argv)<2:
  print("Usage: find_corrupt.py <list of directories>")

corruptlist = []

dirlist = sys.argv[1:]
cflag = False
for directory in dirlist:
  checklist = [fi for fi in os.listdir(directory) if ".post.stdout" in fi]
  for fi in checklist:
    with open("{}/{}".format(directory,fi),'r') as openf:
      for line in openf:
        if "detected corruption" in line:
          corruptlist.append(directory)
          cflag = True
          break
    if cflag:
      cflag = False
      break

for corrupt in corruptlist:
  print(corrupt)
