#!/usr/bin/env python

'Tools for continuum analysis (a la Pippard) of bicycle wheels.'

import numpy as np
from bikewheelcalc import BicycleWheel


def calc_buckling_tension(wheel, approx=None, N=20):
    'Find minimum critical tension within first N modes.'

    def calc_Tc_mode_quad(n):
        'Calculate critical tension for nth mode using full quadratic form.'

        CT = GJ + EIw*n**2/R**2

        if 'y_s' in wheel.rim.sec_params:
            y0 = wheel.rim.sec_params['y_c'] -\
                wheel.rim.sec_params['y_s']
        else:
            y0 = 0.0

        kT = float(ns) / (2*np.pi*R*ls)

        A = -2*kT*n**2*ns*pi*rx**2 + (n**4*ns**2*rx**2)/R**2 -\
            2*kT*n**2*ns*pi*ry**2 + (n**4*ns**2*ry**2)/R**2 +\
            (n**4*ns**2*rx**2*ry**2)/R**4 -\
            (n**2*ns**2*y0)/R + 2*kT*ns*pi*R*y0 - (n**2*ns**2*ry**2*y0)/R**3 +\
            (2*n**4*ns**2*ry**2*y0)/R**3 - 2*kT*n**2*ns*pi*y0**2 +\
            (n**4*ns**2*ry**2*y0**2)/R**4

        B = -2*k_bb*n**2*ns*pi + 4*EI*kT*pi**2 + 4*CT*kT*n**2*pi**2 -\
            (2*EI*n**2*ns*pi)/R**2 - (2*CT*n**4*ns*pi)/R**2 +\
            4*kT*k_bb*pi**2*R**2 - 2*k_uu*n**2*ns*pi*rx**2 -\
            (2*CT*n**4*ns*pi*rx**2)/R**4 - (2*EI*n**6*ns*pi*rx**2)/R**4 -\
            2*k_uu*n**2*ns*pi*ry**2 - (2*EI*n**2*ns*pi*ry**2)/R**4 +\
            (4*EI*n**4*ns*pi*ry**2)/R**4 - (2*EI*n**6*ns*pi*ry**2)/R**4 -\
            (2*k_bb*n**2*ns*pi*ry**2)/R**2 - (4*k_ub*n**2*ns*pi*ry**2)/R +\
            4*k_ub*n**2*ns*pi*y0 + (2*CT*n**2*ns*pi*y0)/R**3 -\
            (2*EI*n**4*ns*pi*y0)/R**3 - (4*CT*n**4*ns*pi*y0)/R**3 +\
            2*k_uu*ns*pi*R*y0 - 2*k_uu*n**2*ns*pi*y0**2 -\
            (2*CT*n**4*ns*pi*y0**2)/R**4 - (2*EI*n**6*ns*pi*y0**2)/R**4

        C = -(-((2*EI*n**2*pi)/R**2) - (2*CT*n**2*pi)/R**2 + 2*k_ub*pi*R)**2 +\
            ((2*CT*n**2*pi)/R**3 + (2*EI*n**4*pi)/R**3 + 2*k_uu*pi*R) *\
            ((2*EI*pi)/R + (2*CT*n**2*pi)/R + 2*k_bb*pi*R)

        # Solve for the smaller root
        T_c = (-B - np.sqrt(B**2 - 4*A*C))/(2*A)

        return T_c

    def calc_Tc_mode_lin(n):
        'Calculate critical tension for nth mode using linear approximation.'

        CT = GJ + EIw*n**2/R**2
        mu = CT/EI

        A = 1 + mu*n**2 + k_bb*R**2/EI
        B = n**4 + mu*n**2
        C = 2*n**2*(1 + mu) - k_ub*R**3/EI
        D = mu*n**2*(n**2 - 1)**2

        f_T = n**2 / (n**2 - R/ls)

        T_c = 2*pi*EI/(ns*R**2*n**2*A) * f_T *\
            (A*k_uu*R**4/EI + B*k_bb*R**2/EI + C*k_ub*R**3/EI + D)

        return T_c

    k_s = calc_continuum_stiff(wheel, tension=0.0)
    k_uu = k_s[0, 0]
    k_ub = k_s[0, 3]
    k_bb = k_s[3, 3]

    # shortcuts
    pi = np.pi
    ns = len(wheel.spokes)
    R = wheel.rim.radius
    ls = wheel.spokes[0].length
    EI = wheel.rim.young_mod * wheel.rim.I22
    EIw = wheel.rim.young_mod * wheel.rim.Iw
    GJ = wheel.rim.shear_mod * wheel.rim.I11

    rx = np.sqrt(wheel.rim.I22 / wheel.rim.area)
    ry = np.sqrt(wheel.rim.I33 / wheel.rim.area)

    if approx == 'linear':

        T_cn = [calc_Tc_mode_lin(n) for n in range(2, N+1)]
        T_c = min(T_cn)
        n_c = T_cn.index(T_c) + 2

    elif approx == 'small_mu':

        T_c = 11.875 * GJ/(ns*R**2) * np.power(k_uu*R**4/GJ, 2.0/3.0)
        n_c = np.power(k_uu*R**4/(2*GJ), 1.0/6.0)

    else:

        T_cn = [calc_Tc_mode_quad(n) for n in range(2, N+1)]
        T_c = min(T_cn)
        n_c = T_cn.index(T_c) + 2

    return T_c, n_c


def calc_continuum_stiff(wheel, tension=None):
    'Calculate smeared-spoke stiffness matrix.'

    if tension is not None:
        wheel.apply_tension(tension)

    k_bar = np.zeros((4, 4))

    for s in wheel.spokes:
        k_bar = k_bar + s.calc_k()/(2*np.pi*wheel.rim.radius)

    return k_bar


def calc_Pn_lat(wheel):
    'Lateral Pippard number (ratio of length scale to spoke spacing).'

    k_sp = calc_continuum_stiff(wheel)
    k_uu = k_sp[0, 0]

    n_spokes = len(wheel.spokes)
    EI = wheel.rim.young_mod * wheel.rim.I22
    GJ = wheel.rim.shear_mod * wheel.rim.I11
    cc = (GJ/EI)/(GJ/EI + 1)

    Pn_lat = n_spokes/(8*np.pi*wheel.rim.radius) *\
        np.power(4*EI * cc / k_uu, 0.25)

    return Pn_lat


def calc_Pn_rad(wheel):
    'Radial Pippard number (ratio of length scale to spoke spacing).'

    k_sp = calc_continuum_stiff(wheel)
    k_vv = k_sp[1, 1]

    n_spokes = len(wheel.spokes)
    EI = wheel.rim.young_mod * wheel.rim.I33

    Pn_rad = n_spokes/(2*np.pi*wheel.rim.radius) *\
        np.power(4*EI / k_vv, 0.25)

    return Pn_rad


def calc_lambda_lat(wheel):
    'Calculate lambda = k_uu*R^4/EI_lat'

    k_sp = calc_continuum_stiff(wheel)
    k_uu = k_sp[0, 0]

    return k_uu*wheel.rim.radius**4 / (wheel.rim.young_mod * wheel.rim.I22)


def calc_lambda_rad(wheel):
    'Calculate lambda = k_vv*R^4/EI_rad'

    k_sp = calc_continuum_stiff(wheel)
    k_vv = k_sp[1, 1]

    return k_vv*wheel.rim.radius**4 / (wheel.rim.young_mod * wheel.rim.I33)


def print_continuum_stats(wheel):
    'Print summary information about the wheel.'

    print('lambda (lat) :', calc_lambda_lat(wheel))
    print('lambda (rad) :', calc_lambda_rad(wheel))
    print('R/le (lat)   :', np.power(calc_lambda_lat(wheel), 0.25))
    print('R/le (rad)   :', np.power(calc_lambda_rad(wheel), 0.25))
    print('Pn_lat       :', calc_Pn_lat(wheel))
    print('Pn_rad       :', calc_Pn_rad(wheel))


def mode_stiff(wheel, n, tension=0.0):
    'Calculate stiffness for the nth mode.'

    k_s = calc_continuum_stiff(wheel, tension)
    k_uu = k_s[0, 0]
    k_ub = k_s[0, 3]
    k_bb = k_s[3, 3]

    # shortcuts
    pi = np.pi
    ns = len(wheel.spokes)
    R = wheel.rim.radius
    # l = wheel.spokes[0].length
    EI = wheel.rim.young_mod * wheel.rim.I22
    EIw = wheel.rim.young_mod * wheel.rim.Iw
    GJ = wheel.rim.shear_mod * wheel.rim.I11

    rx = np.sqrt(wheel.rim.I22 / wheel.rim.area)
    ry = np.sqrt(wheel.rim.I33 / wheel.rim.area)

    CT = GJ + EIw*n**2/R**2

    # Shear center coordinate
    if 'y_s' in wheel.rim.sec_params:
        y0 = wheel.rim.sec_params['y_c'] -\
            wheel.rim.sec_params['y_s']
    else:
        y0 = 0.0

    # Nr = ns*tension / (2*pi)
    Nr = np.sum([s.tension*s.n[1] for s in wheel.spokes]) / (2*pi)

    if n == 0:
        U_uu = 2*pi*R*k_uu
        U_ub = 2*pi*R*k_ub
        U_bb = 2*pi*EI/R + 2*pi*R*k_bb + 2*pi*Nr*y0
    else:  # n > 0
        U_uu = pi*EI*n**4/R**3 + pi*CT*n**2/R**3 + pi*R*k_uu \
            - pi*Nr*n**2/R - pi*Nr*n**2*ry**2/R**3

        U_ub = -pi*EI*n**2/R**2 - pi*CT*n**2/R**2 + pi*R*k_ub \
            - pi*Nr*n**2*ry**2/R**2 - pi*Nr*n**2*y0/R

        U_bb = pi*EI/R + pi*CT*n**2/R + pi*R*k_bb\
            + pi*Nr*y0 - pi*Nr*n**2*(rx**2 + ry**2 + y0**2)/R

    # Solve linear system
    K = np.zeros((2, 2))
    K[0, 0] = U_uu
    K[0, 1] = U_ub
    K[1, 0] = U_ub
    K[1, 1] = U_bb

    x = np.linalg.solve(K, np.array([1, 0]))

    # Displacement stiffness
    Kn_u = 1.0 / x[0]

    # Rotation stiffness
    if x[1] == 0.0:
        Kn_p = float('inf')
    else:
        Kn_p = 1.0 / x[1]

    return Kn_u, Kn_p


def lateral_stiffness(wheel, N=20, tension=0.0):
    'Calculate lateral stiffness for a point load in N/m.'

    Fn = np.zeros(N)  # flexibility of nth mode

    for n in range(len(Fn)):
        Fn[n] = 1.0 / mode_stiff(wheel, n, tension)[0]

    K_lateral = 1.0 / sum(Fn)

    return K_lateral


def lateral_stiffness_phi(wheel, N=20, tension=0.0):
    'Calculate the lateral/rotation stiffness P/phi for a point load.'

    Fn = np.zeros(N)  # flexibility of nth mode

    for n in range(len(Fn)):
        Fn[n] = 1.0 / mode_stiff(wheel, n, tension)[1]

    K_lateral_phi = 1.0 / sum(Fn)

    return K_lateral_phi
