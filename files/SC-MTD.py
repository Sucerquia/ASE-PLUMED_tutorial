from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.nvtberendsen import NVTBerendsen
from gpaw import MixerDif, FermiDirac, GPAW
from ase.calculators.plumed import Plumed
from ase.io import read
from ase import Atoms
from ase import units


pos = read("isomer1.xyz").get_positions()

i = open("plumed.dat", "r").read().splitlines()

T = 100
timestep = 5 * units.fs
taut = 50 * units.fs

atoms = Atoms('Ag5', pos)

MaxwellBoltzmannDistribution(atoms, temperature_K=T)

p = atoms.get_momenta()
psum = p.sum(axis=0) / float(len(p))
p = p - psum
atoms.set_momenta(p)

a = 16

atoms.set_cell([a, a, a])
atoms.set_pbc(True)
atoms.center()

atoms.calc = Plumed(GPAW(h=0.2,
                         mode='lcao',
                         basis='pvalence.dz',
                         xc='PBE',
                         spinpol=True,
                         nbands=-4,
                         occupations=FermiDirac(0.05),
                         parallel=dict(augment_grids=True),
                         mixer=MixerDif(beta=0.25, nmaxold=3, weight=50.0),
                         symmetry='off'),
                    input=i,
                    timestep=timestep,
                    atoms=atoms,
                    kT=units.kB*T)

dyn = NVTBerendsen(atoms, timestep, temperature_K=T, taut=taut, fixcm=False,
                   trajectory='trajectory.traj')

dyn.run(50001)
