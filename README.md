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
- `get_AMQ`
- `solve_Q_e`
- `solve_lakes`

### Discharge Potential from Head
Function named as `Phi_from_fi(fi, k, H)` in `functions.py`

#### Calculations
if <img src="/tex/e6d34c6a5a789e1e4a0bbb3139790504.svg?invert_in_darkmode&sanitize=true" align=middle width=52.398261299999994pt height=22.831056599999986pt/>
<p align="center"><img src="/tex/ab1e28cfe267702e0eeab32d0e2898ae.svg?invert_in_darkmode&sanitize=true" align=middle width=71.3766603pt height=32.990165999999995pt/></p>
else
<p align="center"><img src="/tex/a5affb30940f04eaab5396ac815816c3.svg?invert_in_darkmode&sanitize=true" align=middle width=147.96249435pt height=32.990165999999995pt/></p>

### Head from Discharge Potential
Function named as `fi_from_Phi(Phi, k, H)` in `functions.py`

#### Calculations
if <img src="/tex/c53b171674e90f2c5fbc1efd997a59c5.svg?invert_in_darkmode&sanitize=true" align=middle width=91.01435969999999pt height=27.77565449999998pt/>
<p align="center"><img src="/tex/5e4825e9fb5ed6bc3c4912794e2da345.svg?invert_in_darkmode&sanitize=true" align=middle width=72.1871106pt height=39.452455349999994pt/></p>
else
<p align="center"><img src="/tex/a2fc40d114326e9cf1ef660ebfdc0983.svg?invert_in_darkmode&sanitize=true" align=middle width=107.59562055pt height=37.24318125pt/></p>

### z of <img src="/tex/c91091e68f0e0113ff161179172813ac.svg?invert_in_darkmode&sanitize=true" align=middle width=10.28535419999999pt height=14.15524440000002pt/>
Function named as `z_of_chi(chi, nu, z1, z2)` in `functions.py`

#### Calculations
<p align="center"><img src="/tex/4a263826c44f470d2b8d03be31904eb2.svg?invert_in_darkmode&sanitize=true" align=middle width=121.97223389999999pt height=39.452455349999994pt/></p>
<p align="center"><img src="/tex/8c14bf13c4f436b5b6eb697ba69806ee.svg?invert_in_darkmode&sanitize=true" align=middle width=200.76780899999997pt height=32.990165999999995pt/></p>

### <img src="/tex/c91091e68f0e0113ff161179172813ac.svg?invert_in_darkmode&sanitize=true" align=middle width=10.28535419999999pt height=14.15524440000002pt/> of z
Function named as `chi_of_z(z, nu, z1, z2)` in `functions.py`

#### Calculations
<p align="center"><img src="/tex/4037f5dcb8ef72acd613fd5d9570e1e1.svg?invert_in_darkmode&sanitize=true" align=middle width=123.09491535000001pt height=35.45589465pt/></p>
<p align="center"><img src="/tex/e6327bf7a0513fccdebe5d1af604716c.svg?invert_in_darkmode&sanitize=true" align=middle width=195.4562181pt height=18.80246775pt/></p>

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
Generates a matrix of each elemets far field condition <img src="/tex/983483b08b02816f1eff5947365b4627.svg?invert_in_darkmode&sanitize=true" align=middle width=35.11055129999999pt height=22.465723500000017pt/>

<p align="center"><img src="/tex/2a175680823357bce20a679f228981fd.svg?invert_in_darkmode&sanitize=true" align=middle width=178.0539024pt height=39.452455349999994pt/></p>

### KN Matrix and Q Solver
Function named as `solve_Q_e(AM, Phi0, Phi_lake, N, C, z_ref, W, nw, zw, rw, M, nu, z1, z2, a, m, chi_far)` in `functions.py`

#### Calculations
Generates a vector of each elemets difference in discharge potential without the far field condition <img src="/tex/1c2e35b96f5e17293aa7a7488a480325.svg?invert_in_darkmode&sanitize=true" align=middle width=120.03340634999999pt height=25.680359100000008pt/>

<p align="center"><img src="/tex/c38ee8abebb10ebd13fa97b1ef001cbd.svg?invert_in_darkmode&sanitize=true" align=middle width=86.16824204999999pt height=39.452455349999994pt/></p>

Then the discharges for the far field condition are calculated as

<p align="center"><img src="/tex/fe50763364acdc65f896e7e6ae9da3fa.svg?invert_in_darkmode&sanitize=true" align=middle width=114.41051490000001pt height=17.82653235pt/></p>

### Coefficient Solver
Function named as `solve_lakes(Phi_lake, N, Phi0, z_ref, W, nw, zw, rw, M, nu, z1, z2, m, chi_far)` in `functions.py`

#### Calculations
AN iterative solver for the coefficients <img src="/tex/6512cbd0d448700a036bf3a691c37acc.svg?invert_in_darkmode&sanitize=true" align=middle width=16.81517804999999pt height=14.15524440000002pt/> and <img src="/tex/1afcdb0f704394b16fe85fb40c45ca7a.svg?invert_in_darkmode&sanitize=true" align=middle width=12.99542474999999pt height=22.465723500000017pt/>. It solves firstly for all <img src="/tex/1afcdb0f704394b16fe85fb40c45ca7a.svg?invert_in_darkmode&sanitize=true" align=middle width=12.99542474999999pt height=22.465723500000017pt/> and then for each coefficient <img src="/tex/6512cbd0d448700a036bf3a691c37acc.svg?invert_in_darkmode&sanitize=true" align=middle width=16.81517804999999pt height=14.15524440000002pt/> excluding the element who’s coefficients it is solving for. The loop solving for <img src="/tex/1afcdb0f704394b16fe85fb40c45ca7a.svg?invert_in_darkmode&sanitize=true" align=middle width=12.99542474999999pt height=22.465723500000017pt/> then <img src="/tex/6512cbd0d448700a036bf3a691c37acc.svg?invert_in_darkmode&sanitize=true" align=middle width=16.81517804999999pt height=14.15524440000002pt/> is repeated until the maximum difference between two solves are smaller than a given number, typically <img src="/tex/ec4c25dce67266d40679a50bf9d75e70.svg?invert_in_darkmode&sanitize=true" align=middle width=33.26498669999999pt height=26.76175259999998pt/>.

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
