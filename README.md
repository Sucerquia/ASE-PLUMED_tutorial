# ASE for  Metadynamics Simulations

This tutorial shows how to use the [plumed](https://wiki.fysik.dtu.dk/ase/ase/calculators/plumed.html)
calculator of [ase](https://wiki.fysik.dtu.dk/ase/) for implementing 
Well-Tempered Metadynamics and computing Collective Variables from a 
trajectory.

Plumed actually allows several actions besides of what we show here. For further 
description of plumed details, visit [plumed web page](http://www.plumed.org/doc).


The use of other tools of Plumed with ASE is absolutely analogous to what 
it is explained here.

[[_TOC_]]

# Theory

## Collective Variables

In the most of the cases, it is impossible to extract clear information about 
the system of interest by monitoring the coordinates of all atoms directly, 
even more if our system contains many atoms. Instead of that, it
is possible to make the monitoring simpler by defining functions of those 
coordinates that describe the chemical properties that we are interested in. 
Those functions are called Collective Variables (CVs) and allows biasing 
specific degrees of freedom or analyzing how those properties evolve. Plumed 
has numerous CVs already implemented that can be used with ASE. For a 
complete explanation of CVs implemented in Plumed, 
[go to this link](https://www.plumed.org/doc-v2.8/user-doc/html/colvarintro.html).

## Metadynamics

[Metadynamics](10.1038/s42254-020-0153-0) is an enhanced sampling method 
that allows exploring the configuration landscape by adding cumulative bias in 
terms of some CVs. This bias is added each $`\tau`$ time lapse and usually its 
shape is Gaussian. In time t, the accumulated bias is
defined as:

```math
V_{B}({\bf{s}}, t) = \sum_{t'=\tau, 2\tau,...}^{t'<t}W(t') 
                            \hspace{0.1cm}
                            exp\left({-\sum_i\frac{[s_i\hspace{0.1cm} - 
                            \hspace{0.1cm}
                            s_i(t')]^2}{2\sigma_i}}\right)
```

Where **s** is a set of collective variables, $`\sigma_i`$ is the width of the 
Gaussian related with the i-th collective variable, and *W(t')* is the height 
of the Gaussian in time *t'*. In simple metadynamics, *W(t')* is a constant, 
but in Well-Tempered Metadynamics, the height of the Gaussians is lower where 
previous bias was added. This reduction of height of the new Gaussians decreases 
the error and avoids exploration towards high free energy states that are 
thermodynamically irrelevant. The height in time t' for Well-Tempered 
Metadynamics is defined as:

```math
W(t') = W exp\left({-\frac{\beta \hspace{0.1cm} V_B({\bf s}, 
                  \hspace{0.1cm}t')}{\gamma}}\right)
```

Here, *W* is the maximum height of the Gaussians, $`\beta`$ is the inverse of
the thermal energy ($`1/k_BT`$) and $`\gamma`$ is a bias factor greater than
one that regulates how fast the height of the bias decreases. The higher the 
bias factor, the slower is the decreasing of the heights. Note that when 
$`\gamma`$ approaches infinity, this equation becomes constant and simple 
metadynamics is recovered. On contrast, when $`\gamma`$ approaches zero, no bias is
added, which is the case of Molecular Dynamics.

In this way, the force with bias acting over the i-th atom becomes in:

```math
{\bf F^B}_i = {\bf F}_i - \frac{\partial {\bf s}}{\partial {\bf R}_i} 
                   \frac{\partial V_B({\bf s}, t)}{\partial {\bf s}}
```

where $`{\bf F}_i`$ is the force over the atom i due to interactions with 
the rest of atoms and the second term is the additional force due to the added 
bias.

Part of the power of metadynamics is that it can be used for exploring 
conformations and the accumulated bias converges to the free energy surface 
($`F({\bf s})`$). In the case of Well-Tempered Metadynamics:

```math
\lim_{t\rightarrow \infty} V_B ({\bf{s}}, t) = -\frac{(\gamma -1)}
                                               {\gamma} F({\bf s})
```

# Planar 7-Atoms Cluster

Let's consider a simple system formed by seven atoms with Lennard-Jones (LJ) 
interactions in a planar space. This simple model is presented in the 
[Plumed Masterclass 21.2](https://www.plumed.org/doc-v2.7/user-doc/html/masterclass-21-2.html#masterclass-21-2-ex-9).
This LJ cluster has several stable isomers (Figure 1), which can be 
distinguished in a space of the CVs second (SCM) and third (TCM) central 
moments of the distribution of coordinations (red stars in Figure 2).

![cluster](/cluster.png)

The n-th central moment, $`\mu_n`$, of the coordination number of an N-atoms
in the cluster is defined as

```math
{\mu_n} = \frac{1}{N} \sum_{i=1}^{N} \left( {X}_{i} - 
                \left< {X} \right> \right)^n
```

where $`\left< {X} \right>`$ is the mean value of $`X_i`$, which is the
coordination of the i-th atom:

```math
X_i= \sum_{i\ne j}\frac{1-(r_{ij}/d)^8}{1-(r_{ij}/d)^{16}}
```

For this example, d is fixed to 1.5 $`\sigma`$, in LJ units.

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
Molecular dynamics code:

```python
from ase import units
from ase.io import read
from ase.md.langevin import Langevin
from ase.constraints import FixedPlane
from ase.calculators.plumed import Plumed
from ase.calculators.lj import LennardJones

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
coordinates are imported from `isomer.xyz`, and the 
system remains moving around that state; it does not jump to the other 
isomers. This means we do not obtain a complete sampling of possible 
configurations as mentioned before. Then, an alternative to observe transitions
is to use an enhanced sampling method. In this case, we implement Well-Tempered 
Metadynamics.

```
**Warning**  
Note that in the plumed set-up, there is a line with the keyword UNITS,
which is necessary because all parameters in the plumed set-up and output 
files are assumed to be in plumed internal units. Then, this line is 
important to mantain the units of all plumed parameters and outputs in ASE
units. You can ignore this line if you are aware of the units conversion.
```

### Post Processing Analysis

If you have the trajectory of a MD simulation and you want to compute a set of 
CVs of that trajectory, you can reconstruct the plumed files without running 
again all the simulation. As an example, let's use the trajectory created in 
the last code for rewriting the COLVAR file as follows:

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
         "PRINT ARG=c1.* STRIDE=100 FILE=COLVAR_postpro",
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

As you can see, the first column correspond to the time, the second one is the 
second central moment (SCM) and the third column is the third central moment 
(TCM). When we plot this trajectory in the space of this CVs (that is, the 
second and third columns) we obtain this result:

![](/metadynamics/MD.png)

Where we marked the points that corresponds with the different isomers. Note 
that the system remains confined in the same stable state. That means, for this
case, MD is not enough for exploring all possible configurations and obtaining 
a statistical study of the possible configurations of the system in the 
simulation time scale. Then, an alternative is to use an enhanced sampling 
method. In this case, we implement Well-Tempered Metadynamics for 
reconstructing the Free Energy Surface (FES).

## Well-Tempered Metadynamics Simulation

Well-Tempered Metadynamics method is described in the `Theory section`. It 
basically adds external energy for pushing the system to explore different 
conformations. This makes necessary to add a restrain to avoid that the extra 
energy dissolves the atoms in vacuum. This restriction consists in a 
semi-harmonic potential with this form:

```math
V(d_i)=\left\{
          \begin{array}{ll}
              100 (d_i - 2)^2 & \text{if }d_i>2 \\
              0               & \text{otherwise}
          \end{array}
          \right.
```

Where $`d_i`$ is the distance of each atom to the center of mass. Note that this 
potential does not do anything whereas the distance between the atom and the 
center of mass is lower than 2 (in LJ dimensionless reduced units), but if it 
is greater (trying to escape), this potential begins to work and send it back 
to be close the other atoms. This is defined with the keyword UPPER_WALLS in 
the plumed set up.

Well-Tempered Metadynamics simulation for this case can be run using this 
script:

```python
from ase import units
from ase.io import read
from ase.md.langevin import Langevin
from ase.constraints import FixedPlane
from ase.calculators.plumed import Plumed
from ase.calculators.lj import LennardJones

timestep = 0.005
ps = 1000 * units.fs

setup = [f"UNITS LENGTH=A TIME={1/ps} ENERGY={units.mol/units.kJ}",
         "COM ATOMS=1-7 LABEL=com"]
         "DISTANCE ATOMS=1,com LABEL=d1",
         "UPPER_WALLS ARG=d1 AT=2.0 KAPPA=100.",
         "DISTANCE ATOMS=2,com LABEL=d2",
         "UPPER_WALLS ARG=d2 AT=2.0 KAPPA=100.",
         "DISTANCE ATOMS=3,com LABEL=d3",
         "UPPER_WALLS ARG=d3 AT=2.0 KAPPA=100.",
         "DISTANCE ATOMS=4,com LABEL=d4",
         "UPPER_WALLS ARG=d4 AT=2.0 KAPPA=100.",
         "DISTANCE ATOMS=5,com LABEL=d5",
         "UPPER_WALLS ARG=d5 AT=2.0 KAPPA=100.",
         "DISTANCE ATOMS=6,com LABEL=d6",
         "UPPER_WALLS ARG=d6 AT=2.0 KAPPA=100.",
         "DISTANCE ATOMS=7,com LABEL=d7",
         "UPPER_WALLS ARG=d7 AT=2.0 KAPPA=100.",
         "c1: COORDINATIONNUMBER SPECIES=1-7 MOMENTS=2-3" +
         " SWITCH={RATIONAL R_0=1.5 NN=8 MM=16}",
         "METAD ARG=c1.* HEIGHT=0.05 PACE=500 " +
         "SIGMA=0.1,0.1 GRID_MIN=-1.5,-1.5 GRID_MAX=2.5,2.5" +
         " GRID_BIN=500,500 BIASFACTOR=5 FILE=HILLS"]

atoms = read('isomer.xyz')
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
according to equation `\ref{hills}`. Then, it is necessary to define the 
kT argument of the calculator in ASE units. SIGMA and PACE are the 
standard deviation of the Gaussians and the deposition interval in terms of 
number of steps ($`\tau`$ in equation `\ref{bias}`). HEIGHT and 
BIASFACTOR are the maximum height of the Gaussians (W) and the $`\gamma`$ factor 
of the equation `\ref{hills}`, respectively.

In this case, the Lennard-Jones calculator computes the forces between atoms,
namely, $`{\bf F}_i`$ forces in equation `\ref{bias-force}`. 
Likewise, you could use your preferred ASE calculator.

When one runs a metadynamics simulation, Plumed generates a file called HILLS 
that contains the information of the deposited Gaussians. You can reconstruct 
the free energy by yourself or can use the plumed tool 
[sum_hills](https://www.plumed.org/doc-v2.7/user-doc/html/sum_hills.html). 
The simplest way of using it is:

```
$ plumed sum_hills --hills HILLS
```
After this, Plumed creates a fes.dat file with the FES reconstructed. When the 
FES of this example is plotted, it yields:

![](/metadynamics/fes.png)

## Restart

Suppose you realized it was not enough added bias when it finalized. Then, you 
have to restart your simulation during some steps more. For doing so, you have 
to configure the atoms object in the last state of previous simulation and to 
fix the value of steps in the plumed calculator. Taking the last code as an 
example, this means you would have to change the definition of the object atoms 
as follows:

```python
from ase.io import read
last_configuration = read('MTD.traj')
atoms.set_positions(last_configuration.get_positions())
atoms.set_momenta(last_configuration.get_momenta())
```

and the definition of the calculator becomes in:

```python
atoms.calc = Plumed( ... , restart=True)
atoms.calc.istep = 10000
```

Alternatively, you can initialize your calculator using the next script

```python
from ase.calculators.plumed import restart_from_trajectory

...

atoms.calc = restart_from_trajectory(prev_traj='MTD.traj',
                                     prev_steps = 10000,
                                     ... )
```
