.. _metadynamics:

=================================
ASE for  Metadynamics Simulations
=================================

This tutorial shows how to use the :mod:`~ase.calculators.plumed` calculator
for computing Collective Variables from a Molecular Dynamics (MD) trajectory.
Besides, you will find an example of a Well-Tempered Metadynamics simulation 
using this calculator.

Plumed actually allows several actions like implementing enhanced sampling 
methods (besides metadynamics), to compute other different Collective variables
that are not shown here, among other possibilities. For larger description of 
plumed actions and details, visit `this page <http://www.plumed.org/doc>`_.
The use of other tools of Plumed with ASE is absolutely analogous to what it is
explained here.

Theory
------

Collective Variables
====================

In the most of the cases, it is impossible to extract clear information about 
the system of interest in monitoring the coordinates of all atoms directly, 
even more if our system contains a large quantity of atoms. Instead of that, it
is possible to make the monitoring simpler by defining functions of those 
coordinates that describe the chemical properties that we are interested in. 
Those functions are called Collective Variables (CVs) and allows biasing 
specific degrees of freedom or to analyze how those properties evolve.

An example of CV is the radius of gyration defined by:


.. math::
   \begin{align}
     R= \left(\frac{\sum_i^N |r_i - r_{CM}|^2}{N}\right)^{1/2}
   \end{align}

Where `r_{CM}` is the center of mass of the atoms, *N* is the number of atoms
and `r_i` is the position of atom i. From this definition, it is clear that the
radius of gyration contains information about how disperse is the system in 
respect to the center of mass. That from just one scalar!

Plumed has numerous CVs already implemented that can be used with ASE. For a 
complete explanation of CVs implemented in Plumed, go to 
`here <https://www.plumed.org/doc-v2.7/user-doc/html/colvarintro.html>`_.

Metadynamics
============


:doi:`Metadynamics <10.1038/s42254-020-0153-0>` is an enhanced sampling method 
that allows exploring the configuration landscape by adding cumulative bias in 
terms of some CVs. This bias is added each `\tau` time lapse and usually its 
shape is Gaussian. In time t, that accumulated bias is
defined as:


.. math:: 
   \begin{align}
       V_{B}({\bf{s}}, t) = \sum_{t'=\tau, 2\tau,...}^{t'<t}W(t') 
                            \hspace{0.1cm}
                            exp\left({-\sum_i\frac{[s_i\hspace{0.1cm} - 
                            \hspace{0.1cm}
                            s_i(t')]^2}{2\sigma_i}}\right) 
                            \label{bias} \tag{1}
   \end{align}


Where **s** is a set of collective variables, `\sigma_i` is the width of the 
Gaussian related with the i-th collective variable, and *W(t')* is the height 
of the Gaussian in time *t'*. In simple metadynamics, *W(t)* is a constant, 
but in Well-Tempered Metadynamics, the height of the Gaussians is lower where 
previous bias was added. This reduction of height of the new Gaussians reduces 
the error and avoids exploration towards high free energy states that are 
thermodynamically irrelevant. Then, the height in time t' for Well-Tempered 
Metadynamics is defined as:

.. math:: 
   \begin{align}
       W(t') = W exp\left({-\frac{\beta \hspace{0.1cm} V_B({\bf s}, 
                  \hspace{0.1cm}t')}{\gamma}}\right) \label{hills} \tag{2}
   \end{align}

Here, *W* is the maximum height of the Gaussians, `\beta` is the inverse of
the thermal energy (:math:`1/k_BT`) and `\gamma` is a bias factor greater than
zero that regulates how fast the height of the bias decreases. The higher the 
bias factor, the slower is the decreasing of the heights. Note that when 
`\gamma` approaches infinity, this equation becomes constant and simple 
metadynamics is recovered. On contrast, when `\gamma` approaches 0, no bias is
added, namely, the case of Molecular Dynamics.

In this way, the force acting over the i-th atom becomes in:

.. math::
   \begin{align}
     {\bf F^B}_i = {\bf F}_i - \frac{\partial {\bf s}}{\partial {\bf R}_i} 
                   \frac{\partial V_B({\bf s}, t)}{\partial {\bf s}}
                   \label{bias-force} \tag{3}
   \end{align}

where :math:`{\bf F}_i` is the force over the atom i due to interactions with 
the rest of atoms and the other term is the additional force due to the added 
bias.

Part of the power of metadynamics is that it can be used for exploring 
conformations and the accumulated bias converges to the free energy surface. 
In the case of Well-Tempered Metadynamics, it is:

.. math::
   \begin{align}
      \lim_{t\rightarrow \infty} V_B ({\bf{s}}, t) = -\frac{(\gamma -1)}
                                                      {\gamma} F({\bf s})
   \end{align}


Free Energy Reconstruction
--------------------------

Let's consider a simple system formed by seven atoms with Lennard-Jones (LJ) 
interactions in a planar space. This simple model is presented in the 
`Plumed Masterclass 21.2 <https://www.plumed.org/doc-v2.7/user-doc/html/masterclass-21-2.html#masterclass-21-2-ex-9>`_.
The LJ cluster has several stable isomers (shown below), which can be 
distinguished in a space of the CVs second and third central moments of the 
distribution of coordination numbers.

.. image:: cluster.png
   :width: 800
   :align: center

The nth central moment `\mu_n` of an N-atoms cluster is defined as

.. math::
   \begin{equation}
      {\mu_n} = \frac{1}{N} \sum_{i=1}^{N} \left( {X}_{i} - 
                \left< {X} \right> \right)^n
   \end{equation}

where `X_i` is the coordination number of the i-th atom:

.. math::
   \begin{equation}
      X_i= \sum_{i\ne j}\frac{1-(r_{ij}/d)^8}{1-(r_{ij}/d)^{16}}
   \end{equation}

For this example, d is fixed to `1.5\sigma`, in LJ units.

Molecular Dynamics Simulation
=============================

For showing that it is necessary to use an enhanced sampling method,
let's start with a Langevin simulation without bias. In LJ dimensionless 
reduced units, the parameters of the simulation are  `k_\text{B}T=0.1`, 
friction coefficient fixed equal to 1. 

It is supposed that the system explores all the space of configurations 
due to thermal fluctuations. However, we can see that the system remains in the 
same state, even when we simulate an evolution of the system for a long time 
lapse. That is because a complete exploration of the configuration space could 
take much more time than the possible to simulate. This can be shown running 
the next code (this could take some minutes running):

.. literalinclude:: MD.py

This simulation is started with the configuration of minimum energy, whose 
coordinates are imported from :download:`isomer.xyz`. Note that in the plumed 
set-up, it is added a line with the keyword UNITS. This is necessary because 
all parameters in the plumed set-up are assumed to be in plumed internal units. 
Then, this line is important to remain the units of the plumed parameters in 
the same units as in ASE. You can ignore this line but be aware of the units 
changes.

Post Processing Analysis
========================

If you have the trajectory of a MD simulation and you want to compute a set of 
CVs of that trajectory, you can reconstruct the plumed files without running 
again all the simulation. As an example, let's use the trajectory created in 
the last code for rewriting the COLVAR file as follows:

.. literalinclude:: postpro.py


This code, as well as the previous one, generates a file called COLVAR with 
the value of the CVs. All plumed files begin with a head that describes the 
fields that it contains. In this case, it generates this result::

   $ head -n 2 COLVAR
   #! FIELDS time c1.moment-2 c1.moment-3
    0.000000 0.757954 1.335796

As you can see, the first column correspond to the time, the second one is the 
second central moment (SCM) and the third column is the third central moment 
(TCM). When we plot this trajectory in the space of this CVs (that is, the 
second and third columns) we obtain this result:

.. image:: MD.png
   :width: 400
   :align: center

Where we marked the points that corresponds with the different isomers. Note 
that the system remains confined in the same stable state. That means, for this
case, MD is not enough for exploring all possible configurations and obtaining 
a statistical study of the possible configurations of the system in the 
simulation time scale. Then, an alternative is to use an enhanced sampling 
method. In this case, we implement Well-Tempered Metadynamics for 
reconstructing the Free Energy Surface (FES).


Well-Tempered Metadynamics Simulation
=====================================


This method is described in the theory part. It basically adds external energy 
for pushing the system to explore different conformations. This extra energy 
can generate that the atoms dissolve in vacuum. For that reason, it is 
necessary to add a restrain to avoid it. This restriction consists in add a 
semi-harmonic potential with this form:

.. math::
   \begin{equation}
       V(d_i)=\left\{
          \begin{array}{ll}
              100 (d_i - 2)^2 & \text{if }d_i>2 \\
              0               & \text{otherwise}
          \end{array}
          \right.
   \end{equation}

Where `d_i` is the distance of each atom to the center of mass. Note that this 
potential does not do anything whereas the distance between the atom and the 
center of mass is lower than 2 `\sigma`, but if it is greater (trying to escape),
this potential begins to work and send it back to be close the other atoms.

Then, the code for running this Well-Tempered Metadynamics in ASE is this one:

.. literalinclude:: MTD.py

Note that Well-Tempered Metadynamics requires the value of the temperature 
according to equation :math:`\ref{hills}`. Then, it is necessary to define the 
kT argument of the calculator in the desired units. SIGMA and PACE are the 
initial  standard deviation and the deposition interval in terms of number of 
steps (`\tau` of the equation :math:`\ref{bias}`). HEIGHT and BIASFACTOR are the 
maximum height of the Gaussians (W) and the `\gamma` factor of the equation 
:math:`\ref{hills}`, respectively. 

In this case, the Lennard-Jones calculator computes the forces between atoms, 
namely, :math:`{\bf F}_i` forces in equation :math:`\ref{bias-force}`. 
Likewise, you could use your preferred calculator. Plumed adds the bias forces.

When one runs a metadynamics simulation, Plumed generates a file called HILLS 
that contains the information of the deposited Gaussians. You can reconstruct 
the free energy by yourself or can use the plumed tool 
`sum_hills <https://www.plumed.org/doc-v2.7/user-doc/html/sum_hills.html>`_. 
The simplest way of using it is::

    $ plumed sum_hills --hills HILLS

After this, Plumed creates a fes.dat file with the FES reconstructed. When the 
FES of this example is plotted, it yields:

.. image:: fes.png
   :width: 400
   :align: center

Restart
=======

Suppose you realized it was not enough added bias when it finalized. Then, you 
have to restart your simulation during some steps more. For doing so, you have 
to configure the atoms object in the last state of previous simulation and to 
fix the value of steps in the plumed calculator. Taking the last code as an 
example, this means you would have to change the definition of the object atoms 
as follows::

    from ase.io import read
    last_configuration = read('MTD.traj')
    atoms.set_positions(last_configuration.get_positions())
    atoms.set_momenta(last_configuration.get_momenta())

This allows to 
and the definition of the calculator becomes in::
    
    atoms.calc = Plumed( ... , restart=True)
    atoms.calc.istep = 500000

Alternatively, you can use the next function:

.. autoclass:: ase.calculators.plumed.restart_from_trajectory


