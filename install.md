# Install

In order to use PLUMED with ASE, you need to install the python version of plumed and a version of ASE >= 3.23.0. For specific information about the installation of each code, you might need to check the formal instructions of [PLUMED]( https://www.plumed.org/doc-v2.8/user-doc/html/_installation.html#installingpython ) and [ASE] ( https://wiki.fysik.dtu.dk/ase/install.html#installation-from-source ). However, for a quick installation that guaranties to run this tutorial, you can execute the next lines in your linux terminal if you have conda already installed:

```
# Create an environment
conda create -name plumed-ase
conda activate plumed-ase

# Install py-plumed
conda install -c conda-forge py-plumed -y

# Install ASE
pip install ase
```

##### [Theory: CVs and MTD &rarr;](theory.md)
