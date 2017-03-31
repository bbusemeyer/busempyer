import process_record as pr
import json
import multiprocessing as mp
import sys
from cat_jsons import cat_jsons

# Parallel only worth it for RDM analysis.
def analyze_records(reclist,parallel=True):
  if not parallel:
    arecs = [output_analysis(rec) for rec in reclist]
  else: # Parallelize.
    with mp.Pool(8) as pool:
      arecs = pool.map(output_analysis,reclist)
  return arecs

def output_analysis(recfn):
  print("Processing %s ... "%recfn,end="")
  arec = pr.process_record(json.load(open(recfn,'r')))
  arecfn = recfn.replace("record.json","data.json")
  with open(arecfn,'w') as outf: 
    json.dump(arec,outf)
  print("done.")
  return arecfn

if __name__=="__main__":
  if len(sys.argv)-1 < 1:
    raise AssertionError("""
    Not enough arguements.
    Usage: python analyze_records.py <list of autogen record jsons>
    Make sure form is *record.json.
    """)
  outfn = input("Output file: ")
  reclist = sys.argv[1:]
  areclist = analyze_records(reclist)
  with open(outfn,'w') as outf:
    cat_jsons(areclist,outf)
