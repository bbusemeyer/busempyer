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
  results['ecutwfc'] = float(qeinput.find('basis').find('ecutwfc').text)
  results['ecutrho'] = float(qeinput.find('basis').find('ecutrho').text)
  results['spin_polarized'] = qeinput.find('spin').find('lsda').text == 'true'

  try: results['degauss'] = float(qeinput.find('bands').find('smearing').attrib['degauss'])
  except AttributeError: results['degauss'] = None

  # Structural properties (final if there was a relaxation).
  results['structure'] = read_positions(inpxml)
  results['cell'] = numpy.array([avec.text.split() for avec in qeoutput.find('atomic_structure').find('cell')],dtype=float)

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

  # The convergence indicator seems pretty confused in QE, seems better to compare accuracy to input setting.
  results['converged'] = 'true'==qeoutput.find('convergence_info').find('scf_conv').find('convergence_achieved').text
  results['scf_error'] = float(qeoutput.find('convergence_info').find('scf_conv').find('scf_error').text)

  # Only here if there is smearing and metallic, or something.
  if qeoutput.find('band_structure').find('fermi_energy') is not None:
    results['ef'] = float(qeoutput.find('band_structure').find('fermi_energy').text)
  else:
    results['ef'] = 0.0
  bandstructure = qeoutput.find('band_structure')
  results['nelec'] = int(float(bandstructure.find('nelec').text))
  kptroots  = bandstructure.findall('ks_energies')
  results['kpoints'] = [numpy.array(kptroot.find('k_point').text.split(),dtype=float) for kptroot in kptroots]
  results['bands']  = numpy.array([kptroot.find('eigenvalues').text.split() for kptroot in kptroots],dtype=float) - results['ef']
  results['timing'] = float(root.find('timing_info').find('total').find('wall').text)

  return results

def read_positions(inpxml):
  ''' Read QE XML file into python dict.'''
  tree = parse(inpxml)
  positions = tree.getroot().find('output').find('atomic_structure').find('atomic_positions').findall('atom')
  return {
      'positions':numpy.array([pos.text.split(' ') for pos in positions],dtype=float),
      'species':[pos.attrib['name'] for pos in positions]
    }

def read_kpts(inpxml):
  ''' Read QE XML file into python dict.'''
  tree = parse(inpxml)
  root = tree.getroot()
  kpts = root.find('output').find('band_structure').findall('ks_energies')
  results = [(numpy.array(kpt.find('k_point').text.split(),dtype=float),float(kpt.find('k_point').attrib['weight'])) for kpt in kpts]

  return results

if __name__=='__main__': main()
