# Restart Note

Suppose you realize it was not enough added bias when your simulation finalized. Therefore, you 
have to restart your simulation during some more steps. For doing so, you have 
to configure the atoms object in the last state of previous simulation and to 
fix the value of steps in the plumed calculator. Taking the last code as an 
example, this means you would have to change the initialization of the object atoms 
as follows:

```python
from ase.io import read
last_configuration = read('MTD.traj')
atoms.set_positions(last_configuration.get_positions())
atoms.set_momenta(last_configuration.get_momenta())
```

and the definition of the calculator becomes in

```python
atoms.calc = Plumed( ... , restart=True)
atoms.calc.istep = 10000
```

where the three points must be replaced by the other arguments of the calculator.
Alternatively, you can initialize your calculator using the next script

```python
from ase.calculators.plumed import restart_from_trajectory

...

atoms.calc = restart_from_trajectory(prev_traj='MTD.traj',
                                     prev_steps = 10000,
                                     ... )
```

##### [&larr; Biased simulation: Well Tempered Metadynamics](MTD.md)
##### [Ab-initio: Small Silver Cluster &rarr;](SC.md)