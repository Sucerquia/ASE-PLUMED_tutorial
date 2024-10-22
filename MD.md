## Molecular Dynamics Simulation

For showing that is necessary to use an enhanced sampling method,
let's start with a Langevin simulation without bias. In LJ dimensionless 
reduced units (assuming $`\epsilon`$ = 1 eV, $`\sigma`$ = 1 $`\textrm Å`$ and 
m = 1 a.m.u), the parameters of the simulation are  $`k_\text{B}T=0.1`$, 
friction coefficient fixed equal to 1 and a time step of 0.005.

It is supposed that the system should explore all the space of configurations 
due to thermal fluctuations. However, we can see that the system remains in the 
same state, even when we sample for a long time lapse. That is because a 
complete exploration of the configuration 
space could take more time than the possible to simulate. Figure 2 -blue 
dots- shows the trajectory obtained from the following unbiased 
Molecular dynamics [code](https://gitlab.com/Sucerquia/ase-plumed_tutorial/-/blob/main/files/MD.py):

```python
from ase.calculators.lj import LennardJones
from ase.calculators.plumed import Plumed
from ase.constraints import FixedPlane
from ase.md.langevin import Langevin
from ase.io import read
from ase import units


timestep = 0.005
ps = 1000 * units.fs

setup = [f"UNITS LENGTH=A TIME={1/ps} ENERGY={units.mol/units.kJ}",
         "c1: COORDINATIONNUMBER SPECIES=1-7 MOMENTS=2-3" +
         " SWITCH={RATIONAL R_0=1.5 NN=8 MM=16}",
         "PRINT ARG=c1.* STRIDE=100 FILE=COLVAR",
         "FLUSH STRIDE=1000"]

atoms = read('isomer.xyz')
cons = [FixedPlane(i, [0, 0, 1]) for i in range(7)]
atoms.set_constraint(cons)
atoms.set_masses([1, 1, 1, 1, 1, 1, 1])

atoms.calc = Plumed(calc=LennardJones(rc=2.5, r0=3.),
                    input=setup,
                    timestep=timestep,
                    atoms=atoms,
                    kT=0.1)

dyn = Langevin(atoms, timestep, temperature_K=0.1/units.kB, friction=1,
               fixcm=False, trajectory='UnbiasMD.xyz')

dyn.run(100000)

```

This simulation starts from the configuration of minimum energy, whose 
coordinates are imported from [`isomer.xyz`](https://gitlab.com/Sucerquia/ase-plumed_tutorial/-/blob/main/files/isomer.xyz), and the 
system remains moving around that state; it does not jump to the other 
isomers. This means we do not obtain a complete sampling of possible 
configurations as mentioned before. Then, an alternative to observe transitions
is to use an enhanced sampling method. In this case, we implement Well-Tempered 
Metadynamics.
 
 **NOTE**  
If you want to use the typical `plumed.dat` file to set up the plumed actions, you can replace the `setup` variable assignament in the previous code for:

```python
setup = open("plumed.dat", "r").read().splitlines()
```

| **WARNING** |
| ---         |
| Note that in the plumed set-up, there is a line with the keyword UNITS, which is necessary because all parameters in the plumed set-up and output files are assumed to be in plumed internal units. Then, this line is important to mantain the units of all plumed parameters and outputs in ASE units. You can ignore this line if you are aware of the units  conversion.   |

### Post Processing Analysis

If you have the trajectory of an MD simulation and you want to compute a set of 
CVs of that trajectory, you can reconstruct the plumed files without running 
again all the simulation. As an example, let's use the trajectory created in 
the last code for rewriting the COLVAR file with the next [script](https://gitlab.com/Sucerquia/ase-plumed_tutorial/-/blob/main/files/postpro.py):

```python
from ase.calculators.idealgas import IdealGas
from ase.calculators.plumed import Plumed
from ase.io import read
from ase import units


traj = read('UnbiasMD.traj', index=':')

atoms = traj[0]

timestep = 0.005
ps = 1000 * units.fs
setup = [f"UNITS LENGTH=A TIME={1/ps} ENERGY={units.mol/units.kJ}",
         "c1: COORDINATIONNUMBER SPECIES=1-7 MOMENTS=2-3" +
         " SWITCH={RATIONAL R_0=1.5 NN=8 MM=16}",
         "PRINT ARG=c1.* STRIDE=100 FILE=COLVAR",
         "FLUSH STRIDE=1000"]

calc = Plumed(calc=IdealGas(),
              input=setup,
              timestep=timestep,
              atoms=atoms,
              kT=0.1)

calc.write_plumed_files(traj)
```

This code, as well as the previous one, generates a file called COLVAR with 
the value of the CVs. All plumed files begin with a head that describes the 
fields that it contains. In this case, it generates this result:

```
$ head -n 2 COLVAR
#! FIELDS time c1.moment-2 c1.moment-3
0.000000 0.757954 1.335796
```

As you can see, the first column corresponds to the time, the second one is the 
second central moment (SCM) and the third column is the third central moment 
(TCM). When we plot this trajectory in the space of this CVs (that is, the 
second and third columns) we obtain this result:

<img src="/files/MD.png"  width="400">

**Figure 2.** Unbiased MD trajectory (blue dots) in the space of the collective
variables second and third central moment. Orange stars represent the location of
the local minima isomers of the LJ cluster in this space.

Note that the system remains confined in the same stable state. That means, for this
case, MD is not enough for exploring all possible configurations and obtaining 
a statistical study of the possible configurations of the system -at least in the 
simulated time scale-. Therefore, an alternative is to use an enhanced sampling 
method. In this case, we implement Well-Tempered Metadynamics for 
reconstructing the Free Energy Surface (FES).

##### [&larr; Toy model: Planar 7-Atoms Cluster](defsystem.md)
##### [Biased simulation: Well Tempered Metadynamics &rarr;](MTD.md)
