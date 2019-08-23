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
<p align="center"><img src="/tex/c035769e91775a7eb428b3f25bcb08c4.svg?invert_in_darkmode&sanitize=true" align=middle width=552.1264732499999pt height=67.1688996pt/></p>\underset{uni}{\Omega} = -Wz<p align="center"><img src="/tex/5279c0c08f641f1dba429df6ac010f6c.svg?invert_in_darkmode&sanitize=true" align=middle width=539.3614314pt height=36.164383199999996pt/></p>\underset{well}{\Omega} = \text{NaN}<p align="center"><img src="/tex/cd4a6854dd3bb4ebd2f7f6ef37c0cde1.svg?invert_in_darkmode&sanitize=true" align=middle width=28.242099599999996pt height=11.4155283pt/></p>\underset{well}{\Omega} = \frac{Q}{2\pi}\log(z - z_w)<p align="center"><img src="/tex/5d9b3c813737a840bd281062109b1b3c.svg?invert_in_darkmode&sanitize=true" align=middle width=593.1937506pt height=35.251144499999995pt/></p>\underset{lake}{\Omega} = \text{NaN}<p align="center"><img src="/tex/d56c8ea1f620eae353d8b75eaf46a0c8.svg?invert_in_darkmode&sanitize=true" align=middle width=25.662169499999997pt height=11.4155283pt/></p>\underset{lake}{\Omega} = \sum_2^n \frac{a_n}{\chi^{-n}} + \frac{Q}{2\pi}\log\left(\frac{\chi}{\chi_{far}}\right)$$

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
