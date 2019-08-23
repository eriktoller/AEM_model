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

### Omega Total
Function named as `Omega_total(z, C, W, nw, zw, rw, Q, M, nu, z1, z2, a, m, chi_far, M_not)` in `functions.py`

#### Calculations
<p align="center"><img src="/tex/b6d00e42d6bdb5dc24efdd8e0f29751e.svg?invert_in_darkmode&sanitize=true" align=middle width=176.89787115pt height=21.96341895pt/></p>

### Omega Uniform Flow
Function named as `Omega_uni(z, W)` in `functions.py`

#### Calculations
<p align="center"><img src="/tex/44c030f582ac4cddc27dc379893ea030.svg?invert_in_darkmode&sanitize=true" align=middle width=81.42799335pt height=21.586699199999998pt/></p>

### Omega Well
Function named as `Omega_well(z, zw, rw, Q)` in `functions.py`

#### Calculations
if <img src="/tex/074746e0e5a47c594ef363aed64574c6.svg?invert_in_darkmode&sanitize=true" align=middle width=95.03039699999998pt height=24.65753399999998pt/>
<p align="center"><img src="/tex/167c2937d0714a2eed55f7d25c1e526a.svg?invert_in_darkmode&sanitize=true" align=middle width=79.29816179999999pt height=21.96341895pt/></p>  
else  
<p align="center"><img src="/tex/6ffa722d1184652ed2d6d9324f5e8679.svg?invert_in_darkmode&sanitize=true" align=middle width=152.04851804999998pt height=33.62942055pt/></p>

### Omega Lake
Function named as `Omega_lake(chi, a, Q, chi_far, m)` in `functions.py`

#### Calculations
if <img src="/tex/71acaaa86f32d6fbb06d7e20fc052bb8.svg?invert_in_darkmode&sanitize=true" align=middle width=50.707529399999984pt height=21.18721440000001pt/>  

<p align="center"><img src="/tex/7619f2ee4811aaafcecaba16e469a855.svg?invert_in_darkmode&sanitize=true" align=middle width=79.651506pt height=21.96341895pt/></p>  

else  

<p align="center"><img src="/tex/a738b4421b79a4f6e97224a1adff56ef.svg?invert_in_darkmode&sanitize=true" align=middle width=236.18018159999997pt height=44.69878215pt/></p>

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
