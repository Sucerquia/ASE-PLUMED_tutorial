from ase.calculators.idealgas import IdealGas
from ase.calculators.plumed import Plumed
from ase.io import read
from ase import units


traj = read('UnbiasMD.traj', index=':')

atoms = traj[0]

timestep = 0.005
ps = 1000 * units.fs
setup = open("plumedLJ.dat", "r").read().splitlines()

# IdealGas is a calculator that consider all interactions equal to zero.
calc = Plumed(calc=IdealGas(),
              input=setup,
              timestep=timestep,
              atoms=atoms,
              kT=0.1)

calc.write_plumed_files(traj)