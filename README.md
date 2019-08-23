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
<p align="center"><img src="/tex/80de4694585cb09373e7237104fcd930.svg?invert_in_darkmode&sanitize=true" align=middle width=198.0509223pt height=23.0593242pt/></p>

### Omega Uniform Flow
Function named as `Omega_uni(z, W)` in `functions.py`

#### Calculations
<p align="center"><img src="/tex/f9ed04fc578fd29159c24a1b01900809.svg?invert_in_darkmode&sanitize=true" align=middle width=102.58104449999999pt height=22.68260445pt/></p>

### Omega Well
Function named as `Omega_well(z, zw, rw, Q)` in `functions.py`

#### Calculations
if <img src="/tex/074746e0e5a47c594ef363aed64574c6.svg?invert_in_darkmode&sanitize=true" align=middle width=95.03039699999998pt height=24.65753399999998pt/>
<p align="center"><img src="/tex/8fb0697898d0f7492d4ab22050daf767.svg?invert_in_darkmode&sanitize=true" align=middle width=100.4512146pt height=23.0593242pt/></p>  
else  
<p align="center"><img src="/tex/0d82a5ad0756c30e318f5da5c46a6d75.svg?invert_in_darkmode&sanitize=true" align=middle width=173.2015692pt height=33.62942055pt/></p>

### Omega Lake
Function named as `Omega_lake(chi, a, Q, chi_far, m)` in `functions.py`

#### Calculations
if <img src="/tex/71acaaa86f32d6fbb06d7e20fc052bb8.svg?invert_in_darkmode&sanitize=true" align=middle width=50.707529399999984pt height=21.18721440000001pt/>  

<p align="center"><img src="/tex/96192b498b7b5f331600d4550b52dd3d.svg?invert_in_darkmode&sanitize=true" align=middle width=100.8045588pt height=23.0593242pt/></p>  

else  

<p align="center"><img src="/tex/a5d53d2556b03658023741863dacd162.svg?invert_in_darkmode&sanitize=true" align=middle width=257.33323440000004pt height=44.69878215pt/></p>

### Cauchy Integral
Function named as `Cauchy_integral(N, C, W, nw, zw, rw, Q, M, nu, z1, z2, a, m, chi_far, M_not)` in `functions.py`

#### Calculations
<p align="center"><img src="/tex/5a63bbd852232dc1c7eab122366452c3.svg?invert_in_darkmode&sanitize=true" align=middle width=191.23433175pt height=48.18280005pt/></p>

### AM Matrix
Function named as `get_AMQ(N, z_ref, nw, zw, rw, Q, M, nu, z1, z2, m, chi_far)` in `functions.py`

#### Calculations
<p align="center"><img src="/tex/62d3770aae0f91364540648ee6c200ab.svg?invert_in_darkmode&sanitize=true" align=middle width=103.02125625pt height=39.452455349999994pt/></p>

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
