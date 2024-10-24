
## Well-Tempered Metadynamics Simulation

Well-Tempered Metadynamics method is described in the [`Theory section`](theory.md#metadynamics). It 
basically adds external energy for pushing the system to explore different 
conformations. This makes necessary to add a restrain to avoid that the extra 
energy dissolves the atoms in vacuum. This restriction consists in a 
semi-harmonic potential with the form

```math
V(d_i)=\left\{
          \begin{array}{ll}
              k (d_i - r_w)^2 & \text{if }d_i>2 \\
              0               & \text{otherwise}
          \end{array}
          \right. ~,
```

where $`k`$ is a constraint, in this case being 1000; $`r_w`$ is the location of
the constraint, that we choose to be 2 (in LJ dimensionless reduced units);
$`d_i`$ is the distance of each atom to the center of mass. Note that this 
potential does not do anything whereas the distance between the atom and the 
center of mass is lower than $`r_w`$, but if it 
is greater (trying to escape), this potential begins to work and send it back 
to be close the other atoms. This is defined with the keyword UPPER_WALLS in 
the plumed set up.

You need to create a plumed file with the information of the constraints (also called walls) and the coordinates you want to bias. This file should look like this ([plumedMTD-LJ.dat](https://github.com/Sucerquia/ASE-PLUMED_tutorial/blob/master/files/plumedMTD-LJ.dat)):

```plumed
UNITS LENGTH=A TIME=0.0101805 ENERGY=96.4853329
COM ATOMS=1-7 LABEL=com
DISTANCE ATOMS=1,com LABEL=d1
UPPER_WALLS ARG=d1 AT=2.0 KAPPA=100.
DISTANCE ATOMS=2,com LABEL=d2
UPPER_WALLS ARG=d2 AT=2.0 KAPPA=100.
DISTANCE ATOMS=3,com LABEL=d3
UPPER_WALLS ARG=d3 AT=2.0 KAPPA=100.
DISTANCE ATOMS=4,com LABEL=d4
UPPER_WALLS ARG=d4 AT=2.0 KAPPA=100.
DISTANCE ATOMS=5,com LABEL=d5
UPPER_WALLS ARG=d5 AT=2.0 KAPPA=100.
DISTANCE ATOMS=6,com LABEL=d6
UPPER_WALLS ARG=d6 AT=2.0 KAPPA=100.
DISTANCE ATOMS=7,com LABEL=d7
UPPER_WALLS ARG=d7 AT=2.0 KAPPA=100.
c1: COORDINATIONNUMBER SPECIES=1-7 MOMENTS=2-3 SWITCH={RATIONAL R_0=1.5 NN=8 MM=16}
METAD ARG=c1.* HEIGHT=0.05 PACE=500 SIGMA=0.1,0.1 GRID_MIN=-1.5,-1.5 GRID_MAX=2.5,2.5 GRID_BIN=500,500 BIASFACTOR=5 FILE=HILLS
```

Well-Tempered Metadynamics simulation for this case can be run using the code 
[`MTD.py`](https://gitlab.com/Sucerquia/ase-plumed_tutorial/-/blob/main/files/MTD.py):

```python
from ase.calculators.lj import LennardJones
from ase.calculators.plumed import Plumed
from ase.constraints import FixedPlane
from ase.md.langevin import Langevin
from ase.io import read
from ase import units

timestep = 0.005
ps = 1000 * units.fs

setup = open("plumedMTD-LJ.dat", "r").read().splitlines()

atoms = read('isomerLJ.xyz')
cons = [FixedPlane(i, [0, 0, 1]) for i in range(7)]
atoms.set_constraint(cons)
atoms.set_masses([1, 1, 1, 1, 1, 1, 1])

atoms.calc = Plumed(calc=LennardJones(rc=2.5, r0=3.0),
                    input=setup,
                    timestep=timestep,
                    atoms=atoms,
                    kT=0.1)

dyn = Langevin(atoms, timestep, temperature_K=0.1/units.kB, friction=1,
               fixcm=False, trajectory='MTD.traj')

dyn.run(500000)
```

Note that Well-Tempered Metadynamics requires the value of the temperature 
according to [equation (2)](theory.md#hills) . Then, it is necessary to define the 
kT argument of the calculator. SIGMA and PACE are the 
standard deviation of the Gaussians and the deposition interval in terms of 
number of steps ($`\tau`$ in [equation (1)](theory.md#bias)). HEIGHT and 
BIASFACTOR are the maximum height of the Gaussians (W) and the $`\gamma`$ factor 
of [equation (2)](theory.md#hills), respectively.

In this case, the Lennard-Jones calculator computes the forces between atoms,
namely, $`{\bf F}_i`$ forces in [equation (3)](theory.md#force). 
Likewise, you could use your preferred ASE calculator instead of the LJ calculator used here.

When one runs a metadynamics simulation, Plumed generates a file called HILLS 
that contains the information of the deposited Gaussians. You can reconstruct 
the free energy by yourself or can use the plumed tool 
[sum_hills](https://www.plumed.org/doc-v2.7/user-doc/html/sum_hills.html). 
The simplest way of using it is:

```
$ plumed sum_hills --hills HILLS
```
After this, Plumed creates a fes.dat file with the FES reconstructed. When the 
FES of this example is plotted using [plotter.py](https://gitlab.com/Sucerquia/ase-plumed_tutorial/-/blob/main/files/plotter.py), it yields:

<div align="center">
  <img src="/files/fes.png"  width="400">
</div>

**Figure 3.** Free energy surface (LJ units) in the space of the collective
variables second and third central moment. Orange stars represent the location of
the local minima isomers of the LJ cluster in this space.

Note that theere is bias added around all the states, that means that the system
jumped from one state to the others and gave a complete reconstruction of the
free energy surface.

##### [&larr; Unbiased simulation and Postprocessing](MD.md)
##### [Restart a simulation &rarr;](restart.md)
##### [Ab-initio: Small Silver Cluster &rarr;](SC.md)