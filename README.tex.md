# AEM model
An AEM model program with built-in functions for: head specified boundaries, head specified wells and uniform flow.

**Author: Erik Toller**  
Current version: 0.1 (2019-08-21)  
Master script: `run_aem.py`

## Functions
This groundwater model is constructed out of mutiple complex discharge fucntions with which the solutions are super-positionable. In the current program the following discharge potential functions are included:
- `Omega_total`
- `Omega_uni`
- `Omega_well`
- `Omega_lake`
- `Cauchy_integral`

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

$$\underset{lake}{\Omega}(z) = \sum_2^n \frac{a_n}{\chi^{-n}} + \frac{Q}{2\pi}\log\left(\frac{\chi}{\chi_{far}}\right)$$

### Cauchy Integral
Function named as `Cauchy_integral(N, C, W, nw, zw, rw, Q, M, nu, z1, z2, a, m, chi_far, M_not)` in `functions.py`

#### Calculations
$$a_n = \frac{1}{N}\sum_k^N \underset{other}{\Omega}(z_k)\e^{-in\theta_k}$$

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
