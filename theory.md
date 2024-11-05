# Theory

## Collective Variables

In most cases, it is impossible to extract clear information about 
a system of interest by monitoring the Cartesian coordinates of all atoms directly, 
especially if our system contains many atoms. Instead, it
is possible to monitor the system by defining functions of those 
coordinates that describe the chemical properties that we are interested in. 
These functions are called Collective Variables (CVs) and allow biasing 
specific degrees of freedom or analyzing how those properties evolve. Plumed 
has numerous CVs already implemented that can be used with ASE. For a 
complete explanation of CVs implemented in Plumed, 
[go to this link](https://www.plumed.org/doc-v2.8/user-doc/html/colvarintro.html).

## Metadynamics

[Metadynamics](https://www.nature.com/articles/s42254-020-0153-0) is an enhanced sampling method 
that allows exploring the configuration landscape by adding cumulative bias in 
terms of some CVs. This bias is added each $`\tau`$ time lapse and usually its 
shape is Gaussian. In time t, the accumulated bias is
defined as

<a name="bias"></a>
$$
(1)
$$

$$
V_{B}({\mathbf{s}}, t) = \sum_{t'=\tau, 2\tau,...}^{t' < t}W(\mathbf{s}, t') 
                            \hspace{0.1cm}
                            \exp\left({-\sum_i\frac{[s_i\hspace{0.1cm} - 
                            \hspace{0.1cm}
                            s_i(t')]^2}{2\sigma_i}}\right) ~,
$$

where $`\mathbf{s}`$ is a set of collective variables, $`\sigma_i`$ is the width of the 
Gaussian related to the i-th collective variable, and *W(s, t')* is the height 
of the Gaussian in time *t'*. In simple metadynamics, *W(s, t')* is a constant, 
but in Well-Tempered Metadynamics, the height of the Gaussians is lower where
previous bias was added. This reduction in the height of the new Gaussians decreases 
the error and avoids exploration towards high free energy states that are 
thermodynamically irrelevant. The height in time t' for Well-Tempered 
Metadynamics is

<a name="hills"></a>
$$
(2)
$$

$$
W(\mathbf{s}, t')  = W \exp\left({-\frac{\beta \hspace{0.1cm} V_B({\bf s}, 
                  \hspace{0.1cm}t')}{\gamma}}\right) ~,
$$

with $`W`$ the maximum height of the Gaussians, $`\beta`$ the inverse of
the thermal energy ($`1/k_BT`$) and $`\gamma`$ a bias factor greater than
one that regulates how fast the height of the bias decreases: the higher the 
bias factor, the slower the heights decrease. Note that when 
$`\gamma`$ approaches infinity, this equation becomes constant and simple 
metadynamics is recovered. In contrast, when $`\gamma`$ approaches zero, no bias is
added, which is the case of Molecular Dynamics.

The addition of the bias potential produces an extra force in each atom, such that
the resultant force for the i-th atom is

<a name="force"></a>
$$
(3)
$$

$$
{\bf F^B}_i = {\bf F}_i - \frac{\partial {\bf s}}{\partial {\bf R}_i} 
                   \frac{\partial V_B({\bf s}, t)}{\partial {\bf s}} ~,
$$

where $`{\bf F}_i`$ is the natural unbiased force, $`R_i`$ is the coordinate of the atom, and the second term 
is the additional force due to the added bias.

Part of the power of metadynamics is that it can be used for exploring 
conformations. Moreover, the accumulated bias converges to the free energy surface 
($`F({\bf s})`$),

$$
(4)
$$

$$
\lim_{t\rightarrow \infty} V_B ({\bf{s}}, t) = -\frac{(\gamma -1)}
                                               {\gamma} F({\bf s})~.
$$

##### [&larr; Introduction](NAVIGATION.md)
##### [Toy model: Planar 7-Atoms Cluster &rarr;](defsystem.md)
