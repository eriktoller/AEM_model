# -*- coding: utf-8 -*-
"""
Function files for AEM program

This is a collection of defined functions for an Analytical Element Method
solution program.

Erik Toller
2019-07-27
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt


def Phi_from_fi(fi, k, H):
    if fi < H:
        Phi = 0.5*k*fi**2
    else:
        Phi = k*H*fi - 0.5*k*H**2
    return Phi


def fi_from_Phi(Phi, k, H):
    cond = 0.5*k*H**2
    if Phi < cond:
        fi = np.sqrt(2*Phi/k)
    else:
        fi = (Phi + 0.5*k*H**2) / (k*H)
    return fi


def Omega_uni(z, W):
    Omega = -W * z
    return(Omega)


def Omega_well(z, zw, rw, Q):
    if abs(z - zw) < rw:
        Omega = complex(np.nan, np.nan)
    else:
        Omega = Q / (2 * cmath.pi) * cmath.log(z - zw)
    return(Omega)


def Omega_total(z, C, W, nw, zw, rw, Q,
                M, nu, z1, z2, a, m, chi_far, M_not):
    Omega = C + 0j
    if W != 0:
        Omega += Omega_uni(z, W)
    if nw != 0:
        for ii in range(nw):
            Omega += Omega_well(z, zw[ii], rw[ii], Q[0, M+ii])
    if M > 0:
        for mm in range(M):
            if mm != M_not:
                chi = chi_of_z(z, nu[mm], z1[mm], z2[mm])
                Omega += Omega_lake(chi, a[mm, :],
                                    Q[0, mm], chi_far[mm], m[mm])
    return Omega


def z_of_chi(chi, nu, z1, z2):
    try:
        z1 = z1[0]
        z2 = z2[0]
        nu = nu[0]
    except TypeError:
        z1 = z1
        z2 = z2
        nu = nu
    Z = .5 * (chi / nu + nu / chi)
    z = .5 * (Z * (z2 - z1) + z2 + z1)
    return z


def chi_of_z(z, nu, z1, z2):
    try:
        z1 = z1[0]
        z2 = z2[0]
        nu = nu[0]
    except TypeError:
        z1 = z1
        z2 = z2
        nu = nu
    Z = (2 * z - z1 - z2) / (z2 - z1)
    Zm = Z - complex(1, 0)
    Zp = Z + complex(1, 0)
    if np.imag(Zm) == 0:
        Zm = np.real(Zm)
    elif np.imag(Zp) == 0:
        Zp = np.real(Zp)
    chi = nu * (Z + cmath.sqrt(Zm) * cmath.sqrt(Zp))
    return chi


def Omega_lake(chi, a, Q, chi_far, m):
    chi_con = np.conj(chi)
    if chi * chi_con < 0.999:
        Omega = complex(np.nan, np.nan)
    else:
        if m == 0:
            Omega = complex(0, 0)
        else:
            Omega = a[m]
            for ii in range(m):
                mm = m-(ii+1)
                if mm != 0:
                    Omega = Omega/chi + a[mm]
                else:
                    Omega /= chi
    if Q != 0:
        Omega += Q / (2 * cmath.pi) * cmath.log(chi / abs(chi_far))
    return Omega


def Cauchy_integral(N, C, W, nw, zw, rw, Q,
                    M, nu, z1, z2, a, m, chi_far, M_not):
    d_theta = 2 * cmath.pi / N
    theta0 = 0.5*d_theta
    intgrl = np.zeros([N, m[M_not[0]]+1], dtype=np.complex_)
    b = np.zeros([1, m[M_not[0]]+1], dtype=np.complex_)
    for ii in range(N):
        theta = theta0 + ii * d_theta
        chi = cmath.exp(1j*theta)
        z = z_of_chi(chi, nu[M_not[0]], z1[M_not[0]], z2[M_not[0]])
        Omega = Omega_total(z, C, W, nw, zw, rw, Q,
                            M, nu, z1, z2, a, m, chi_far, M_not[-1])
        for jj in range(m[M_not[0]]+1):
            intgrl[ii, jj] = np.real(Omega) * cmath.exp(-1j * jj * theta)
    for kk in range(m[M_not[0]]+1):
        b[0, kk] = 2 / N * sum(intgrl[:, kk])
    b[0, 0] /= 2
    return b


def solve_Q_e(AM, Phi0, Phi_lake, N, C, z_ref, W, nw, zw, rw,
              M, nu, z1, z2, a, m, chi_far):
    KN = np.zeros([M+nw+1, 1])
    Q0 = np.zeros([1, M+nw])
    C = 0  # <-- This one has been missing from the ellipse funciton
    for ii in range(M):
        aa = Cauchy_integral(N, C, W, nw, zw, rw, Q0,
                             M, nu, z1, z2,
                             a, m, chi_far, [ii, -1])
        KN[ii, 0] = Phi_lake[ii] - np.real(aa[0, 0])
    for jj in range(nw):
        aa = Omega_total(zw[jj]+rw[jj]*1.0001, C, W, nw, zw, rw, Q0,
                         M, nu, z1, z2, a, m, chi_far, -1)
        KN[M+jj, 0] = Phi_lake[M+jj] - np.real(aa)
    KN[M+nw, 0] = Phi0 - np.real(Omega_total(z_ref, 0, W, nw, zw, rw, Q0,
                                             M, nu, z1, z2, a, m, chi_far, -1))
    Q = np.transpose(np.linalg.lstsq(AM, KN)[0])
    return Q


def get_AMQ(N, z_ref, nw, zw, rw, Q, M, nu, z1, z2, m, chi_far):
    AM = np.zeros([M+nw+1, M+nw+1])
    a = np.zeros([M, max(m)+1])
    for ii in range(M):
        for ij in range(M):
            Q *= 0
            Q[0, ij] = 1
            m = [0]*M
            C = 0
            W = 0
            aa = Cauchy_integral(N, C, W, nw, zw, rw, Q,
                                 M, nu, z1, z2,
                                 a, m, chi_far, [ii, -1])
            AM[ii, ij] = np.real(aa[0, 0])
        for ik in range(nw):
            zk = z_of_chi(1, nu[ii], z1[ii], z2[ii])
            AM[ii, M+ik] = np.real(Omega_well(zk, zw[ik], rw[ik], 1))
        AM[ii, M+nw] = 1
    for jj in range(nw):
        for ji in range(M):
            aa = cmath.log(chi_of_z(zw[jj], nu[ji], z1[ji], z2[ji])
                           / abs(chi_far[ji])) / (2*cmath.pi)
            AM[M+jj, ji] = np.real(aa)
        for jk in range(nw):
            AM[M+jj, M+jk] = np.real(Omega_well((zw[jj]+rw[jj]*1.0001),
                                                zw[jk], rw[jk], 1))
        AM[M+jj, M+nw] = 1
    for kk in range(M):
        aa = cmath.log(chi_of_z(z_ref, nu[kk], z1[kk], z2[kk]) /
                       abs(chi_far[kk])) / (2*cmath.pi)
        AM[M+nw, kk] = np.real(aa)
    for ki in range(nw):
        AM[M+nw, M+ki] = np.real(Omega_well(z_ref, zw[ki], rw[ki], 1))
    AM[M+nw, M+nw] = 1
    return AM


def solve_lakes(Phi_lake, N, Phi0, z_ref, W, nw, zw, rw,
                M, nu, z1, z2, m, chi_far):
    error = 1
    NIT = 0
    C = Phi0
    a = np.zeros([M, max(m)+1], dtype=np.complex_)
    a_old = np.zeros([M, max(m)+1], dtype=np.complex_)
    a_solve = np.zeros([M, max(m)+1], dtype=np.complex_)
    Q = np.zeros([1, M+nw+1])
    Q_old = np.zeros([1, M+nw+1])
    AM = get_AMQ(N, z_ref, nw, zw, rw, Q, M, nu, z1, z2, m, chi_far)
    while error > 1e-6 or not NIT < 200:
        Q = solve_Q_e(AM, Phi0, Phi_lake, N, C, z_ref, W, nw, zw, rw,
                      M, nu, z1, z2, a, m, chi_far)
        C = Q[0, -1]
        for ii in range(M):
            a_solve[ii, 0:(m[ii]+1)] = Cauchy_integral(N, C, W, nw, zw, rw, Q,
                                                       M, nu, z1, z2, a, m,
                                                       chi_far, [ii])
            a[ii, :] = - np.conj(a_solve[ii, :])
        error_a = np.max(abs(a-a_old))
        error_Q = np.max(abs(Q-Q_old))
        if error_a > error_Q:
            error = error_a
        else:
            error = error_Q
        NIT += 1
        a_old[:, :] = a[:, :]
        Q_old = Q
        print('')
        print('Iteration: ', NIT)
        print('')
        print('Error')
        print(error)
        print('Error a')
        print(error_a)
        print('Error Q')
        print(error_Q)
        print('')
    print('--- Solver complete ---')
    print('')
    print('     Iterations:', NIT)
    print('  Maximum error:', error)
    print('Maximum error a:', error_a)
    print('Maximum error Q:', error_Q)
    return a, Q, C


def Contour_flow_net(Nx, xmin, xmax, Ny, ymin, ymax, lvs,
                     C, W, nw, zw, rw, Q,
                     M, nu, z1, z2, a, m, chi_far, M_not, N):
    print('')
    print('Plotting flow net')
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    fig, ax = plt.subplots()
    X = np.linspace(xmin, xmax, num=Nx)
    Y = np.linspace(ymin, ymax, num=Ny)
    Z = np.zeros([Nx, Ny], dtype=np.complex_)
    for xx in range(Nx):
        for yy in range(Ny):
            z = complex(X[xx], Y[yy])
            Z[yy, xx] = Omega_total(z, C, W, nw, zw, rw, Q,
                                    M, nu, z1, z2, a, m, chi_far, M_not)
    if M > 0:
        for jj in range(M):
            d_theta = 2 * cmath.pi / N
            z = [0]*(N+1)
            for ii in range(N):
                theta = ii * d_theta
                chi = cmath.exp(1j*theta)
                z[ii] = z_of_chi(chi, nu[jj], z1[jj], z2[jj])
            z[N] = z[0]
            ax.plot(np.real(z), np.imag(z), color='black', linewidth=0.5)
    CS_real = ax.contour(X, Y, np.real(Z), levels=lvs,
                         colors='red', linewidths=0.5)
    CS_imag = ax.contour(X, Y, np.imag(Z), levels=lvs,
                         colors='blue', linewidths=0.5)
    h1, _ = CS_real.legend_elements()
    h2, _ = CS_imag.legend_elements()
    ax.legend([h1[0], h2[0]], ['Equipotential', 'Stream lines'])
    ax.set_title('Flow Net')


def Contour_head(Nx, xmin, xmax, Ny, ymin, ymax, lvs,
                 C, W, nw, zw, rw, Q,
                 M, nu, z1, z2, a, m, chi_far, M_not, N, k, H):
    print('')
    print('Plotting head contour')
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    fig, ax = plt.subplots()
    X = np.linspace(xmin, xmax, num=Nx)
    Y = np.linspace(ymin, ymax, num=Ny)
    Z = np.zeros([Nx, Ny])
    for xx in range(Nx):
        for yy in range(Ny):
            z = complex(X[xx], Y[yy])
            Omega = Omega_total(z, C, W, nw, zw, rw, Q,
                                M, nu, z1, z2, a, m, chi_far, M_not)
            Z[yy, xx] = fi_from_Phi(np.real(Omega), k, H)
    if M > 0:
        for jj in range(M):
            d_theta = 2 * cmath.pi / N
            z = [0]*(N+1)
            for ii in range(N):
                theta = ii * d_theta
                chi = cmath.exp(1j*theta)
                z[ii] = z_of_chi(chi, nu[jj], z1[jj], z2[jj])
            z[N] = z[0]
            ax.plot(np.real(z), np.imag(z), color='black', linewidth=0.5)
    CS = ax.contour(X, Y, Z, levels=lvs, linewidths=0.5)
    ax.clabel(CS, inline=1, fontsize=10)
    ax.set_title('Hydraulic Head')


def check_heads(C, W, nw, zw, rw, Q, M, nu, z1, z2,
                a, m, chi_far, M_not, N, k, H, fi, z_ref, fi_ref):
    Phi = np.real(Omega_total(z_ref, C, W, nw, zw, rw, Q, M, nu,
                              z1, z2, a, m, chi_far, -1))
    fi1 = fi_from_Phi(Phi, k, H)
    fi_array = [0]*(M+nw)
    d_fi = [0]*(M+nw)
    for ii in range(M):
        z = z_of_chi(1, nu[ii], z1[ii], z2[ii])
        Omega = Omega_total(z, C, W, nw, zw, rw, Q,
                            M, nu, z1, z2, a, m, chi_far, M_not)
        Phi = np.real(Omega)
        fi_M = fi_from_Phi(Phi, k, H)
        fi_array[ii] = fi_M
    for ii in range(nw):
        z = zw[ii]+rw[ii]*1.0001
        Omega = Omega_total(z, C, W, nw, zw, rw, Q,
                            M, nu, z1, z2, a, m, chi_far, M_not)
        Phi = np.real(Omega)
        fi_nw = fi_from_Phi(Phi, k, H)
        fi_array[M+ii] = fi_nw
    for ii in range(M+nw):
        d_fi[ii] = fi[ii] - fi_array[ii]
    print('')
    print('--- Head check ---')
    print('')
    print('Differance in head at:')
    print('Reference is:', fi1 - fi_ref)
    for ii in range(M):
        print('    M =', ii, 'is:', d_fi[ii])
    for ii in range(nw):
        print('   nw =', ii, 'is:', d_fi[M+ii])
    return d_fi
