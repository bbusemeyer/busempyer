
def read_dir_espresso(espinp):
  record = {}

  inpstr = ''
  for line in espinp:
    inpstr += line
  inplines = inpstr.split('\n')

  for lidx,line in enumerate(inplines):
    if '=' in line:
      spl = line.split()
      if len(spl) < 3: 
        # Should be "this = that" not this=that or this =that
        print "Warning: line not formatted right (%s)"%espinp
      try:
        if '.' in spl[2]:
          record[spl[0]] = float(spl[2])
        else:
          record[spl[0]] = int(spl[2])
      except ValueError:
        record[spl[0]] = spl[2]

    if 'K_POINTS' in line:
      record['kpoint'] = map(int,inplines[lidx+1].split())

  return record
