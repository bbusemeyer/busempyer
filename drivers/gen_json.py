#!/usr/bin/python

import json
import sys
import data_processing as dp
from mython import NumpyToListEncoder
from subprocess import check_output
from imp import reload
reload(dp)

# Neat way of calling: 
# find . -name '*_metadata.json' > rootlist
# python gen_json.py $(< rootlist) &> gen_json.out

files = sys.argv[1:]
roots = [f.replace('_metadata.json','') for f in files]
locs = ['/'.join(f.split('/')[:-1]) for f in files]
locs = set(locs)
if len(locs)>1:
  raise AssertionError
else:
  loc = list(locs)[0]
for root in roots:
  data = dp.read_dir_autogen(root,gosling='/home/brian/bin/gosling')
  outfn = loc+"/record.json"
  print("Outputting to %s..."%outfn)
  with open(outfn,'w') as outf:
    json.dump(data,outf,cls=NumpyToListEncoder)
