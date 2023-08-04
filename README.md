# Curlometer Extension
We reformat and extend the Curlometer technique so that it reconstructs magnetic fields using more than 4 spacecraft.

## Methods
Traditionally, the Curlometer technique uses magnetic field measurements from four spacecraft to estimate the current density near the spacecraft configuration. If we assume that bulk plasma velocity is much slower than the speed of light, then we can simplify Ampere's law to be $\nabla \times B = \mu_0 J$. It therefore follows that the current density can be computed from the estimations of $\partial B_m$.

We modify this method to instead estimate the value of $B$ itself at points in space not measured by a spacecraft. This problem is formulated as a first order Taylor Series $\forall i \in$ {1,2,3,4}, $m \in$ {*x,y,z*}

```math
\hat{B}_m^i=B_m+\sum_{k \in \{x,y,z\}} \partial_k B_{m} r_{k}^{i}.
```


## How to Use


# References
This work was originally published in the 2021 Frontiers in Astronomy and Space Sciences article *Magnetic Field Reconstruction for a Realistic Multi-Point, Multi-Scale Spacecraft Observatory* by Broeren et al.
[open access link](https://doi.org/10.3389/fspas.2021.727076)

