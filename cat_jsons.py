def cat_jsons(inplist,outf):
  outf.write('[\n')
  for jfn in inplist[1:-1]:
    outf.write(open(jfn,'r').read())
    outf.write(',\n')
  outf.write(open(inplist[-1],'r').read())
  outf.write('\n]')

if __name__ == "__main__":
  import sys

  if len(sys.argv) < 2:
    raise AssertionError("""
    Incorrect number of arguements.
    Usage: cat_jsons.py <list of jsons>
    """)

  outfn = raw_input('Output file: ')
  inplist = sys.argv[1:-1]
  with open(outfn,'w') as outf:
    cat_jsons(inplist,outf)
