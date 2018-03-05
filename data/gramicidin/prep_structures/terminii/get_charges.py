import parmed as pm
from parmed.charmm import CharmmPsfFile, CharmmCrdFile, CharmmParameterSet


print('Loading CHARMM files...')
params = CharmmParameterSet('toppar/par_all36m_prot.prm', 'toppar/par_all36_lipid.prm')

gram = CharmmPsfFile('../../simulations/themis_data/1jnoetaohh3o2.psf')

for atom in gram.atoms:
    print(atom.name, atom.idx, atom.charge, atom.rmin, atom.epsilon)