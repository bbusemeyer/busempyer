import process_record as pr
import json
import sys

reclist = sys.argv[1:]
if len(reclist) < 2:
  print("Usage: python analyze_records.py <list of autogen record jsons>")
  print("Make sure form is *record.json.")

for recfn in reclist:
  print("Processing %s ... "%recfn,end="")
  arec = pr.process_record(json.load(open(recfn,'r')))
  with open(recfn.replace("record.json","data.json"),'w') as outf: 
    json.dump(arec,outf)
  print("done.")
