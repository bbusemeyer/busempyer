#!/usr/bin/python

import json
import sys
from data_processing import read_dir
from mython import NumpyToListEncoder
from subprocess import check_output

# Neat way of calling: 
# find . -name '*_metadata.json' > rootlist
# python gen_json.py $(< rootlist) &> gen_json.out

files = sys.argv[1:]
roots = [f.replace('_metadata.json','') for f in files]
for root in roots:
  data = read_dir(root,gosling='/home/brian/bin/gosling')
  for k in data.keys():
    print "Dumping {0}_{1} to JSON".format(root,k)
    with open('{0}_{1}_datarecord.json'.format(root,k),'w') as outf:
      json.dump(data[k],outf,cls=NumpyToListEncoder)
