# ASE-PLUMED Tutorial
### A tutorial about the ASE-PLUMED calculator presented in [Sucerquia *et. al.*, JCP, 2022](https://doi.org/10.1063/5.0082332) [[1](#cite1)]
by Daniel Sucerquia, Pilar Cossio and Olga Lopez-Acevedo.

This tutorial shows how to use the [plumed calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/plumed.html)
in [ASE](https://wiki.fysik.dtu.dk/ase/).

The Atomic Simulation Environment, [ASE](https://wiki.fysik.dtu.dk/ase/), is a package that allows setting up,
run and visualize atomistic simulations. It is interfaced with some [other codes](https://wiki.fysik.dtu.dk/ase/#supported-calculators),
which use quantum or classical methods.  We developed a [plumed calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/plumed.html) to connect PLUMED with all the other codes interfaced in ASE.  We note that PLUMED allows several actions in addition of what we show in this tutorial. For further description, visit [plumed web page](http://www.plumed.org/doc). 

This tutorial begins with a brief explanation of the basic ideas about metadynamics. We first use a toy model to show how to compute collective variables on-the-fly in an MD simulation, and also from a previously computed MD trajectory using post-processing. We then implement Well-Tempered Metadynamics to reconstruct the free energy surface of the toy system along a set of collective variables. Finally, we present a more realistic example for a small silver cluster that uses ab-initio metadynamics to estimate the free energy landscape. All the files required to complete this tutorial are publicly available in [https://github.com/Sucerquia/ASE-PLUMED_tutorial/blob/master/files](https://github.com/Sucerquia/ASE-PLUMED_tutorial/blob/master/files).

<a name="cite1">[1]</a> Sucerquia, Daniel, Cristian Parra, Pilar Cossio, and Olga Lopez-Acevedo. “Ab Initio Metadynamics Determination of Temperature-Dependent Free-Energy Landscape in Ultrasmall Silver Clusters.” The Journal of Chemical Physics 156, no. 15 (April 21, 2022): 154301. [https://doi.org/10.1063/5.0082332](https://doi.org/10.1063/5.0082332)

| **WARNING** |
| ---         |
| In order to complete this tutorial you need to install [py-plumed](https://www.plumed.org/doc-v2.8/user-doc/html/_installation.html#installingpython) and a version of [ASE](https://gitlab.com/ase/ase) >= 3.23.0. For further details, check [installation instructions](install.md).|

```mermaid
flowchart TB;
  Z[PLUMED syntax] -.-> A[ASE-PLUMED tutorial]
  A -.-> B[Install packages]
  A ==> C[Theory: CVs and MTD]
  A ==> D[Toy model: Planar 7-Atoms Cluster]
  D ==> E[Unbiased simulation and Postprocessing]
  D ==> F[Biased simulation: Well Tempered Metadynamics]
  F -.-> G[Restart a simulation]
  A ==> H[Ab-initio: Small Silver Cluster]

  click Z "ref1" "A previous tutorial that introduces the basics of PLUMED syntax"
  click B "install.md" "Install py-plumed and ASE"
  click C "theory.md" "Theory necessary to follow completely this tutorial"
  click D "defsystem.md" "Toy model showing the usage of PLUMED in ASE"
  click E "MD.md" "Unbiased simulation computing CVs on the fly and by postprocessing"
  click F "MTD.md" "Addition of the bias using metadynamics"
  click G "restart.md" "How to restart a simulation in ASE"
  click H "SC.md" "Ab-initio example"
```

##### [Theory: CVs and MTD &rarr;](theory.md)
