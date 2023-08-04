# Curlometer Extension
We reformat and extend the Curlometer technique so that it reconstructs magnetic fields using more than 4 spacecraft.

## Methods
Traditionally, the Curlometer technique uses magnetic field measurements from four spacecraft to estimate the current density near the spacecraft configuration. If we assume that bulk plasma velocity is much slower than the speed of light, then we can simplify Ampere's law to be $\nabla \times B = \mu_0 J$. It therefore follows that the current density can be computed from the estimations of $\partial B_m$.

### Reformulation
We modify this method to instead estimate the value of $B$ itself at points in space not measured by a spacecraft. This problem is formulated as a first order Taylor Series $\forall i \in$ {1,2,3,4}, $m \in$ {*x,y,z*}

```math
\hat{B}_m^i=B_m+\sum_{k \in \{x,y,z\}} \partial_k B_{m} r_{k}^{i}.
```

In this equation $\hat{B_m}^i$ is the measured $m^{th}$ component of $B$ at the $i^{th}$ spacecraft, $B_m$ is the computed $m^{th}$ component of $B$ at $\xi$, $\partial_k B_{m}$ is the computed derivative of the $m^{th}$ component of $B$ with respect to the $k^{th}$ direction at $\xi$, and $r_k^{i}$ is the relative position of spacecraft $i$ with respect to $\xi$. In other words, if $x_{ik}$ is the $k^{th}$ component of spacecraft $i$'s location, then $r_k^{i} := x_{ik} - \xi_k$.

This equation can be reformatted into an equation in matrix form with 12 equations and 12 unknowns.

```math
 \begin{bmatrix} \hat{B}_x^{1} \\ \hat{B}_x^{2} \\ \hat{B}_x^{3} \\ \hat{B}_x^{4} \\ \hat{B}_y^{1} \\ \hat{B}_y^{2} \\ \hat{B}_y^{3} \\ \hat{B}_y^{4} \\ \hat{B}_z^{1} \\ \hat{B}_z^{2} \\ \hat{B}_z^{3} \\ \hat{B}_z^{4}  \end{bmatrix}  
=
\begin{bmatrix}
1 & r_x^{1} & r_y^{1} & r_z^{1} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
1 & r_x^{2} & r_y^{2} & r_z^{2} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
1 & r_x^{3} & r_y^{3} & r_z^{3} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
1 & r_x^{4} & r_y^{4} & r_z^{4} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & r_x^{1} & r_y^{1} & r_z^{1} & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & r_x^{2} & r_y^{2} & r_z^{2} & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & r_x^{3} & r_y^{3} & r_z^{3} & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & r_x^{4} & r_y^{4} & r_z^{4} & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & r_x^{1} & r_y^{1} & r_z^{1} \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & r_x^{2} & r_y^{2} & r_z^{2} \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & r_x^{3} & r_y^{3} & r_z^{3} \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & r_x^{4} & r_y^{4} & r_z^{4}
\end{bmatrix}
\begin{bmatrix} B_x \\ \partial_x B_{x} \\ \partial_y B_{x} \\ \partial_z B_{x} \\ B_y \\ \partial_x B_{y} \\ \partial_y B_{y} \\ \partial_z B_{y} \\ B_z \\ \partial_x B_{z} \\ \partial_y B_{z} \\ \partial_z B_{z}  \end{bmatrix} .
```
By solving this linear system for the values in the right hand vector, we get an estimate of all 3 components of $B$ at the reconstructed point.

### Extension
We extend this method to configurations of more than 4 spacecraft by using combinatorics to our advantage. If we have a system of $N$ spacecraft, then it will contain 
```math
\begin{pmatrix} N \\ 4 \end{pmatrix} = \frac{N!}{4! (N-4)!}
```
tetrahedral subsets (i.e groups 4 spacecraft). We preferentially select the tetrahedron of spacecraft that pass both of the following criteria:
1. The subset must have a shape parameter less than 1. We define the shape parameter, $\chi$, in terms of the tetrahedron's elongation, $E$, and planarity, $P$:
```math
\chi = \sqrt{E^2 + P^2} < 1
```
2. The barycenter of the tetrahedron of spacecraft, $r_0$, must fall relatively close to the point in space, $\xi$, that we wish to reconstruct the magnetic field. Nearby is defined as being within one characteristic size, $L$, of the barycenter:
```math
|r_0 - \xi| < L
```

We ignore estimates from tetrahedra that do not pass these criteria. At each reconstructed point, we now estimate the value of $B$ using all tetrahedra that pass these criteria. From this esemble of estimates we take the median value as our final computed magnetic field vector $B$.



## How to Use
We have included a select amount of data from a numerical simulation of plasma turbulence (solar wind like conditions). This data, stored in the folder 'example_data' can be used to reconstruct one example magnetic field around 3 different HelioSwarm (9-spacecraft) configurations. To change which configuration of spacecraft is used, cycle the *hour* parameter between 94, 144, and 205.
```
hour = 205              # hour of HS configuration to select [94, 144, or 205]
```
To change the nearby criteria, change the value of the *L_coeff* parameter. This modifies the equation above so that $|r_0 - \xi| < L_{coeff} L$
```
L_coeff = 1             # radius of reconstruction around each tetrahedron's barycenter
```
To change the shape criteria, change the *chi_thres* parameter. This modifies the equation above so that $\chi = \sqrt{E^2 + P^2} < \chi_{thres}$
```
chi_thres = 1.0         # shape threshold for using a tetrahedron in reconstruction (chi := sqrt{E^2 + P^2})
```



# References
This work was originally published in the 2021 Frontiers in Astronomy and Space Sciences article *Magnetic Field Reconstruction for a Realistic Multi-Point, Multi-Scale Spacecraft Observatory* by Broeren et al.
[open access link](https://doi.org/10.3389/fspas.2021.727076)

