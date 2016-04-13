import data_processing as dp
import json
import sys

reclist = sys.argv[1:]
if len(reclist) < 2:
  print("Usage: python analyze_records.py <list of autogen record jsons>")
  print("Make sure form is <something>.record.json")

for recfn in reclist:
  print("Processing %s ..."%recfn,)
  arec = dp.process_record(json.load(open(recfn,'r')))
  with open(rec.replace("record.json","data.json"),'w') as outf: 
    outf.dump(arec)
  print("done.")
