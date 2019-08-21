# -*- coding: utf-8 -*-
"""
This is a program test for AEM modelling

Erik Toller
2019-07-29
"""

from functions import (Contour_flow_net,
                       Contour_head,
                       solve_lakes,
                       chi_of_z,
                       Omega_total,
                       fi_from_Phi,
                       Phi_from_fi,
                       check_heads)
import matplotlib.pyplot as plt
import numpy as np

plt.close('all')

"Reference point and aquifer properties"
H = 100
k = 1
z_ref = complex(1000, 0)
fi_ref = 4*2

"Well data"
zw = [complex(75, 20), complex(20, -20)]
rw = [0.1, 0.1]
nw = 2
W = complex(.1, 0)*0
fi_nw = [4, 7]

"Lake data"
M = 3
nu = [0.677032961426901, 0.723043208990021, 1, 1, 1]
z1 = [complex(0, 0), complex(-40, -10), complex(-35, 25),
      complex(-35, 25), complex(0, -40)]
z2 = [complex(30, 20), complex(-40, 10), complex(0, -40),
      complex(35, 25), complex(30, 10)]
m = [20]*M
N = 200
fi_M = [8, 12, 9, 9, 9]

"Plot properties"
Nx = 200
Ny = Nx
xmin = -100
xmax = 100
ymin = -100
ymax = 100
lvs = 30

"Pre-calculations"
chi_far = [0]*M
for ii in range(M):
    chi_far[ii] = chi_of_z(z_ref, nu[ii], z1[ii], z2[ii])

fi = fi_M[0:M] + fi_nw[0:nw]
Phi0 = Phi_from_fi(fi_ref, k, H)
Phi_lake = [0]*(M+nw)
for ii in range(M+nw):
    Phi_lake[ii] = Phi_from_fi(fi[ii], k, H)

"Solver for coeficients, discharges and constant"
[a, Q, C] = solve_lakes(Phi_lake, N, Phi0, z_ref, W, nw, zw, rw,
                        M, nu, z1, z2, m, chi_far)

"Checking the heads"
d_fi = check_heads(C, W, nw, zw, rw, Q, M, nu, z1, z2,
                   a, m, chi_far, -1, N, k, H, fi, z_ref, fi_ref)

"Plotters"
Contour_head(Nx, xmin, xmax, Ny, ymin, ymax, lvs,
             C, W, nw, zw, rw, Q,
             M, nu, z1, z2, a, m, chi_far, -1, N, k, H)

Contour_flow_net(Nx, xmin, xmax, Ny, ymin, ymax, lvs,
                 C, W, nw, zw, rw, Q,
                 M, nu, z1, z2, a, m, chi_far, -1, N)
