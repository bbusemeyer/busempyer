#!/usr/bin/python

import json
import sys
import data_processing as dp
import mython as my

# Neat way of calling: 
# find . -name '*.inp' > rootlist
# python gen_json_espresso.py $(< rootlist) &> gen_json.out

files = sys.argv[1:]
roots = [f.replace('.inp','') for f in files]
for root in roots:
  data = dp.read_dir_espresso(root)
  print "Dumping {0} to JSON".format(root)
  with open('{0}_datarecord.json'.format(root),'w') as outf:
    json.dump(data,outf,cls=my.NumpyToListEncoder)
