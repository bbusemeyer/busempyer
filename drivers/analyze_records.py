import process_record as pr
import json
import multiprocessing as mp
import sys
debug = True

reclist = sys.argv[1:]
if len(reclist) < 2:
  print("Usage: python analyze_records.py <list of autogen record jsons>")
  print("Make sure form is *record.json.")

def output_analysis(recfn):
  print("Processing %s ... "%recfn,end="")
  arec = pr.process_record(json.load(open(recfn,'r')))
  with open(recfn.replace("record.json","data.json"),'w') as outf: 
    json.dump(arec,outf)
  print("done.")

if debug:
  for rec in reclist:
    output_analysis(rec)
else: # Parallelize.
  with mp.Pool(8) as pool:
    arecs = pool.map(output_analysis,reclist)
