'''
Homemade cif parser for comparing cif properties. 
Sometimes pymatgen takes liberties with the structure you give it.
This tries to keep things as raw as possbile while still being useful.
'''
def seperr(num):
  try:
    pp = num.find('(')
    ep = num.find(')')
    dp = num.find('.')
    if pp < 0:
      return float(num),0.0
    else:
      return float(num[:pp]), float(num[pp+1:ep])*10**-(pp-dp-1)
  except ValueError:
    return num,0.0

def read_cif(cifstr):
  # Useful keys that have the pattern:
  # _key value
  one_keys=[
      '_cell_length_a','_cell_length_b','_cell_length_c',
      '_cell_angle_alpha','_cell_angle_beta','_cell_angle_gamma',
      '_symmetry_space_group_name_H-M','_symmetry_Int_Tables_number',
      '_cell_measurement_temperature','_cell_measurement_pressure'
    ]

  cifdat={}
  siteflag=False

  for line in ciffn.split('\n'):
    lilist=line.split()

    if len(lilist)==0: continue

    if lilist[0] in one_keys:
      cifdat[lilist[0]],cifdat[lilist[0]+'_err']=seperr(lilist[1])

    if lilist[0]=='#End':
      siteflag=False
    if siteflag:
      cifdat['sites'].append({
          'atom':lilist[0],
          'pos':[seperr(it)[0] for it in lilist[4:7]],
          'pos_err':[seperr(it)[1] for it in lilist[4:7]]
        })
    if lilist[0]=='_atom_site_attached_hydrogens':
      siteflag=True
      cifdat['sites']=[]

  return cifdat

def test_on(ciffn):
  results=read_cif(ciffn)
  #print(results)
  print(results['sites'])

if __name__=='__main__':
  test_on('lowtemp_fese.cif')
