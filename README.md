# AEM model
An AEM model program with built-in functions for: head specified boundaries, head specified wells and uniform flow.

**Author: Erik Toller**  
Current version: 0.1 (2019-08-21)  
Master script: `run_aem.py`

## Functions
This groundwater model is constructed out of mutiple complex discharge fucntions with which the solutions are super-positionable. In the current program the following functions are included:
- `Omega_uni`
- `Omega_well`
- `Omega_lake`

### Omega Uniform Flow
Function named as `Omega_uni(z, W)` in `functions.py`

#### Calculations
<img src="/tex/08c04de91cdce60f5e68be561db33bd5.svg?invert_in_darkmode&sanitize=true" align=middle width=81.42799334999998pt height=22.465723500000017pt/>

### Omega Well
Function named as `Omega_well(z, zw, rw, Q)` in `functions.py`

#### Calculations
if <img src="/tex/074746e0e5a47c594ef363aed64574c6.svg?invert_in_darkmode&sanitize=true" align=middle width=95.03039699999998pt height=24.65753399999998pt/>  

<img src="/tex/edb2f5dceb7043b2f9f302310ab517f5.svg?invert_in_darkmode&sanitize=true" align=middle width=79.29816344999999pt height=22.465723500000017pt/>  

else  

<img src="/tex/0e874dbb95e888192231b27293774abc.svg?invert_in_darkmode&sanitize=true" align=middle width=148.52173049999996pt height=30.392597399999985pt/>

### Omega Lake
Function named as `Omega_lake(chi, a, Q, chi_far, m)` in `functions.py`

#### Calculations
if <img src="/tex/71acaaa86f32d6fbb06d7e20fc052bb8.svg?invert_in_darkmode&sanitize=true" align=middle width=50.707529399999984pt height=21.18721440000001pt/>  

<img src="/tex/3d817b6781989b088a3dde0759b36703.svg?invert_in_darkmode&sanitize=true" align=middle width=79.65150599999998pt height=22.465723500000017pt/>  

else  

$\underset{lake}{\Omega} = \sum_1^m 

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
