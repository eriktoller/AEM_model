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
                       Phi_from_fi,
                       check_heads)
import matplotlib.pyplot as plt

plt.close('all')

"Reference point and aquifer properties"
H = 100
k = 1
z_ref = complex(1000, 0)
fi_ref = 8
W = complex(.1, 0)

"Well data"
zw = [complex(75, 75), complex(20, -20)]
rw = [0.1, 0.1]
fi_nw = [4, 7]
nw = 1

"Lake data"
z1 = [complex(0, 0), complex(-40, -10), complex(-35, 25),
      complex(-35, 25), complex(0, -40)]
z2 = [complex(30, 20), complex(-40, 10), complex(0, -40),
      complex(-5, 0), complex(-5, 0)]
nu = [0.677032961426901, 0.723043208990021, 1, 1, 1]
N = 200
fi_M = [8, 10, 9, 9, 9]
M = 1
if M == 0:
    m = [10]*1
else:
    m = [10]*M

"Impermeable element data"
M_imp = 2
nu_imp = [0.6, 0.5]
z1_imp = [complex(-75, 80), complex(75, -75)]
z2_imp = [complex(-50, 70), complex(75, 0)]
m_imp = [10, 10]


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
[a, Q, C, alpha] = solve_lakes(Phi_lake, N, Phi0, z_ref, W, nw, zw, rw,
                               M, nu, z1, z2, m, chi_far,
                               M_imp, nu_imp, z1_imp, z2_imp, m_imp)

"Checking the heads"
d_fi = check_heads(C, W, nw, zw, rw, Q, M, nu, z1, z2,
                   a, m, chi_far, -1,
                   M_imp, nu_imp, z1_imp, z2_imp, alpha, m_imp, -1,
                   N, k, H, fi, z_ref, fi_ref)

"Plotters"
Contour_head(Nx, xmin, xmax, Ny, ymin, ymax, lvs,
             C, W, nw, zw, rw, Q,
             M, nu, z1, z2, a, m, chi_far, -1,
             M_imp, nu_imp, z1_imp, z2_imp, alpha, m_imp, -1,
             N, k, H)

Contour_flow_net(Nx, xmin, xmax, Ny, ymin, ymax, lvs,
                 C, W, nw, zw, rw, Q,
                 M, nu, z1, z2, a, m, chi_far, -1,
                 M_imp, nu_imp, z1_imp, z2_imp, alpha, m_imp, -1, N)
