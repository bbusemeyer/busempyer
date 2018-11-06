#!/usr/bin/python3
import sys
import os
import subprocess as sub

if len(sys.argv) < 2:
  raise AssertionError("""No input. Usage: 

  python iter_autogen.py <file with list of paths to scripts>

  List in file is delineated by newlines.
  """)
scriptlist = open(sys.argv[1],'r').read().split('\n')

for script in scriptlist:
  if os.path.isfile(script):
    exestr = "python {} &> {} &".format(script,script.replace('.py','.out'))
    print("Running \n {} \n ".format(exestr))
    sub.call( exestr, shell=True )
  else:
    print("Couldn't find {}.".format(script))
#sub.call("tail -f " + \
#    " ".join([ script.replace(".py",".out") for script in scriptlist ]),
#    shell=True)
