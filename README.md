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
<p align="center"><img src="/tex/2099f28d0ed0bb4858b9e66a88f733b8.svg?invert_in_darkmode&sanitize=true" align=middle width=85.81931325pt height=21.83559675pt/></p>

### Omega Well
Function named as `Omega_well(z, zw, rw, Q)` in `functions.py`

#### Calculations
if <img src="/tex/074746e0e5a47c594ef363aed64574c6.svg?invert_in_darkmode&sanitize=true" align=middle width=95.03039699999998pt height=24.65753399999998pt/>
<p align="center"><img src="/tex/a7f43e2a27c1965851c158705d2e2102.svg?invert_in_darkmode&sanitize=true" align=middle width=85.116075pt height=22.0182633pt/></p>  
 <p align="center"><img src="/tex/5610a87adca69423b0f2ac3de333cc70.svg?invert_in_darkmode&sanitize=true" align=middle width=157.8664296pt height=33.62942055pt/></p>

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
