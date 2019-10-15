from xml.etree.ElementTree import parse
import numpy

def main():
  results = read_xml('unrestrict.xml')
  print(results)

def read_xml(inpxml):
  ''' Read QE XML file into python dict.'''
  tree = parse(inpxml)
  root = tree.getroot()
  qeinput = root.find('input')
  qeoutput = root.find('output')

  results = {}
  results['nelec']  = 63*4+5
  results['ecutwfc'] = float(qeinput.find('basis').find('ecutwfc').text)
  results['ecutrho'] = float(qeinput.find('basis').find('ecutrho').text)
  results['spin_polarized'] = qeinput.find('spin').find('lsda').text == 'true'

  try: results['degauss'] = float(qeinput.find('bands').find('smearing').attrib['degauss'])
  except AttributeError: results['degauss'] = None

  # Reading kpoints.
  kibz = qeinput.find('k_points_IBZ').find('monkhorst_pack')
  if kibz is None: results['nkpts'] = (1,1,1)
  else: results['nkpts'] = (int(kibz.attrib['nk1']),int(kibz.attrib['nk2']),int(kibz.attrib['nk3']))

  if qeinput.find('bands').find('smearing') is not None:
    results['degauss'] = float(qeinput.find('bands').find('smearing').attrib['degauss'])

  # Bands (format depends on spin).
  if results['spin_polarized']:
    results['nbndup'] = int(qeoutput.find('band_structure').find('nbnd_up').text)
    results['nbnddn'] = int(qeoutput.find('band_structure').find('nbnd_dw').text)
  else:
    results['nbndup'] = int(qeoutput.find('band_structure').find('nbnd').text)
    results['nbnddn'] = int(qeoutput.find('band_structure').find('nbnd').text)

  results['magtot'] = float(qeoutput.find('magnetization').find('total').text)
  results['magabs'] = float(qeoutput.find('magnetization').find('absolute').text)
  results['energy'] = float(qeoutput.find('total_energy').find('etot').text)

  # Only here if there is smearing and metallic, or something.
  if qeoutput.find('band_structure').find('fermi_energy') is not None:
    results['ef'] = float(qeoutput.find('band_structure').find('fermi_energy').text)
  else:
    results['ef'] = 0.0
  evalroots  = qeoutput.find('band_structure').find('ks_energies').findall('eigenvalues')
  results['bands']  = numpy.array([evalroot.text.split() for evalroot in evalroots],dtype=float) - results['ef']
  results['timing'] = float(root.find('timing_info').find('total').find('wall').text)

  return results

def read_kpts(inpxml):
  ''' Read QE XML file into python dict.'''
  tree = parse(inpxml)
  root = tree.getroot()
  kpts = root.find('output').find('band_structure').findall('ks_energies')
  results = [(numpy.array(kpt.find('k_point').text.split(),dtype=float),float(kpt.find('k_point').attrib['weight'])) for kpt in kpts]

  return results

if __name__=='__main__': main()
