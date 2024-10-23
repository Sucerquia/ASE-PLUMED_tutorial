# Theory

## Collective Variables

In most of cases, it is impossible to extract clear information about 
a system of interest by monitoring the coordinates of all atoms directly, 
even more if our system contains many atoms. Instead of that, it
is possible to make the monitoring simpler by defining functions of those 
coordinates that describe the chemical properties that we are interested in. 
Those functions are called Collective Variables (CVs) and allows biasing 
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
```math
(1)
```
```math
V_{B}({\bf{s}}, t) = \sum_{t'=\tau, 2\tau,...}^{t'<t}W(t') 
                            \hspace{0.1cm}
                            exp\left({-\sum_i\frac{[s_i\hspace{0.1cm} - 
                            \hspace{0.1cm}
                            s_i(t')]^2}{2\sigma_i}}\right) ~,
```

where **s** is a set of collective variables, $`\sigma_i`$ is the width of the 
Gaussian related with the i-th collective variable, and *W(t')* is the height 
of the Gaussian in time *t'*. In simple metadynamics, *W(t')* is a constant, 
but in Well-Tempered Metadynamics, the height of the Gaussians is lower where 
previous bias was added. This reduction of height of the new Gaussians decreases 
the error and avoids exploration towards high free energy states that are 
thermodynamically irrelevant. The height in time t' for Well-Tempered 
Metadynamics is 

<a name="hills"></a>
```math
(2)
```
```math
W(t') = W exp\left({-\frac{\beta \hspace{0.1cm} V_B({\bf s}, 
                  \hspace{0.1cm}t')}{\gamma}}\right) ~,
```

with *W* the maximum height of the Gaussians, $`\beta`$ the inverse of
the thermal energy ($`1/k_BT`$) and $`\gamma`$ a bias factor greater than
one that regulates how fast the height of the bias decreases: the higher the 
bias factor, the slower is the decreasing of the heights. Note that when 
$`\gamma`$ approaches infinity, this equation becomes constant and simple 
metadynamics is recovered. On contrast, when $`\gamma`$ approaches zero, no bias is
added, which is the case of Molecular Dynamics.

The addition of the bias potential produces an extra force in each atom, such that
the resultant force for the i-th atom is

<a name="force"></a>
```math
(3)
```
```math
{\bf F^B}_i = {\bf F}_i - \frac{\partial {\bf s}}{\partial {\bf R}_i} 
                   \frac{\partial V_B({\bf s}, t)}{\partial {\bf s}} ~,
```

where $`{\bf F}_i`$ is the natural unbiased force and the second term 
is the additional force due to the added bias.

Part of the power of metadynamics is that it can be used for exploring 
conformations. Moreover, the accumulated bias converges to the free energy surface 
($`F({\bf s})`$),

```math
(4)
```
```math
\lim_{t\rightarrow \infty} V_B ({\bf{s}}, t) = -\frac{(\gamma -1)}
                                               {\gamma} F({\bf s})~.
```
##### [Toy model: Planar 7-Atoms Cluster &rarr;](defsystem.md)