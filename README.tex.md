# AEM model
An AEM model program with built-in functions for: head specified boundaries, head specified wells and uniform flow.

**Author: Erik Toller**  
Current version: 0.1 (2019-08-21)  
Master script: `run_aem.py`

## Functions
This groundwater model is constructed out of mutiple complex discharge fucntions with which the solutions are super-positionable. In the current program the following discharge potential functions are included:
- `Phi_from_fi`
- `fi_from_Phi`
- `z_of_chi`
- `chi_of_z`
- `Omega_total`
- `Omega_uni`
- `Omega_well`
- `Omega_lake`
- `Cauchy_integral`
- `get_AMQ`
- `solve_Q_e`
- `solve_lakes`

### Discharge Potential from Head
Function named as `Phi_from_fi(fi, k, H)` in `functions.py`

#### Calculations
if $\phi < H$
$$\Phi = \frac{1}{2}k\phi^2$$
else
$$\Phi = kH\phi - \frac{1}{2}kH^2$$

### Head from Discharge Potential
Function named as `fi_from_Phi(Phi, k, H)` in `functions.py`

#### Calculations
if $\Phi < \frac{1}{2}kH^2$
$$\phi = \sqrt{2\frac{\Phi}{k}}$$
else
$$\phi = \frac{\Phi + \frac{1}{2}kH^2}{kH}$$

### z of $\chi$
Function named as `z_of_chi(chi, nu, z1, z2)` in `functions.py`

#### Calculations
$$Z = \frac{1}{2}\left(\frac{\chi}{\nu} + \frac{\nu}{\chi}\right)$$
$$z = \frac{1}{2}(Z(z_2 - z_1) + z_2 + z_1)$$

### $\chi$ of z
Function named as `chi_of_z(z, nu, z1, z2)` in `functions.py`

#### Calculations
$$Z = \frac{2z - z_1 - z_2}{z_2 - z_1}$$
$$\chi = \nu(Z + \sqrt{Z-1}\sqrt{Z+1})$$

### Omega Total
Function named as `Omega_total(z, C, W, nw, zw, rw, Q, M, nu, z1, z2, a, m, chi_far, M_not)` in `functions.py`

#### Calculations
$$\Omega(z) = \underset{uni}{\Omega} + \underset{well}{\Omega} + \underset{lake}{\Omega} + C$$

### Omega Uniform Flow
Function named as `Omega_uni(z, W)` in `functions.py`

#### Calculations
$$\underset{uni}{\Omega}(z) = -Wz$$

### Omega Well
Function named as `Omega_well(z, zw, rw, Q)` in `functions.py`

#### Calculations
if $|z - z_w| < r_w$
$$\underset{well}{\Omega}(z) = \text{NaN}$$  
else  
$$\underset{well}{\Omega}(z) = \frac{Q}{2\pi}\log(z - z_w)$$

### Omega Lake
Function named as `Omega_lake(chi, a, Q, chi_far, m)` in `functions.py`

#### Calculations
if $\chi\bar{\chi} < 1$  

$$\underset{lake}{\Omega}(z) = \text{NaN}$$  

else  

$$\underset{lake}{\Omega}(z) = \sum_2^n \frac{a_n}{\chi^{n}} + \frac{Q}{2\pi}\log\left(\frac{\chi}{\chi_{far}}\right)$$

### Cauchy Integral
Function named as `Cauchy_integral(N, C, W, nw, zw, rw, Q, M, nu, z1, z2, a, m, chi_far, M_not)` in `functions.py`

#### Calculations
$$a_n = \frac{1}{N}\sum_{k=1}^N \underset{other}{\Omega}(z_k)\text{e}^{-in\theta_k}$$

### AM Matrix
Function named as `get_AMQ(N, z_ref, nw, zw, rw, Q, M, nu, z1, z2, m, chi_far)` in `functions.py`

#### Calculations
Generates a matrix of each elemets far field condition $\Lambda_{m,n}$

$$ \text{AM} = \begin{bmatrix} 
\Lambda_{1,1} & \Lambda_{1,n} & 1 \\
\Lambda_{m,1} & \Lambda_{m,n} & 1
\end{bmatrix}$$

### KN Matrix and Q Solver
Function named as `solve_Q_e(AM, Phi0, Phi_lake, N, C, z_ref, W, nw, zw, rw, M, nu, z1, z2, a, m, chi_far)` in `functions.py`

#### Calculations
Generates a vector of each elemets difference in discharge potential without the far field condition $\beta_{m} = \underset{head}{\Phi} - \underset{\text{no Q}}{\Phi}$

$$ \text{KN} = \begin{bmatrix} 
\beta_{1} \\
\beta_{m}
\end{bmatrix}$$

Then the discharges for the far field condition are calculated as

$$Q = \textbf{AM}^{-1}\textbf{KN}$$

### Coefficient Solver
Function named as `solve_lakes(Phi_lake, N, Phi0, z_ref, W, nw, zw, rw, M, nu, z1, z2, m, chi_far)` in `functions.py`

#### Calculations
AN iterative solver for the coefficients $a_n$ and $Q$. It solves firstly for all $Q$ and then for each coefficient $a_n$ excluding the element who’s coefficients it is solving for. The loop solving for $Q$ then $a_n$ is repeated until the maximum difference between two solves are smaller than a given number, typically $10^{-6}$.

## Definition of variables
**Reference point and aquifer properties**  
`H` aquifer thickness  
`k` hydraulic conductivity  
`z_ref` complex coordinate reference point  
`fi_ref` hydraulic head at reference point  
`W` uniform flow (complex)  

**Well data**  
`zw` list of complex coordinates for wells  
`rw` list of radii of wells  
`fi_nw` hydraulic head at wells  
`nw` number of wells  

**Lake data**  
`z1` list of complex coordinates for ellipse locii start  
`z2` list of complex coordinates for ellipse locii end  
`nu` list of &nu; for ellipse  
`m` list of number of coefficient for head specified boundary elements  
`N` number of integration point for Cauchy integral  
`fi_M` hydraulic head at head specified boundary elements  
`M` number of head specified boundary elements  

**Plot properties**  
`Nx` grid plot size in x-direction  
`Ny` grid plot size in y-direction  
`xmin` start value in x-direction  
`xmax` start value in x-direction  
`ymin` start value in y-direction  
`ymax` start value in y-direction  
`lvs` number of contour levels  
