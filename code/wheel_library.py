#!/usr/bin/env python

"""Generate wheels from standard templates or rules"""

import numpy as np
import pandas as pd
from bikewheelcalc import BicycleWheel, Hub, Rim, Spoke


def wheel_from_name(name):
    'Generate a wheel from a standard template'

    df = pd.read_csv('../data/wheel_properties.csv', index_col=0)

    # Find requested wheel
    if name in df:
        p = df[name]

        w = BicycleWheel()
        w.hub = Hub(diam1=float(p['hub_diameter'])/1000.,
                    width1=float(p['hub_width']) / 2000.)
        w.rim = Rim(radius=float(p['rim_radius'])/1000.,
                    area=float(p['rim_area'])/1e6,
                    I11=float(p['rim_GJ']) / (float(p['rim_shear_mod'])*1e9),
                    I22=float(p['rim_EI2']) / (float(p['rim_young_mod'])*1e9),
                    I33=float(p['rim_EI1']) / (float(p['rim_young_mod'])*1e9),
                    Iw=0.0,
                    young_mod=float(p['rim_young_mod'])*1e9,
                    shear_mod=float(p['rim_shear_mod'])*1e9)

        if p['spk_pattern'] == 'radial':
            w.lace_radial(n_spokes=int(p['spk_num']),
                          diameter=float(p['spk_diameter'])/1000.,
                          young_mod=float(p['spk_young_mod'])*1e9,
                          offset=float(p['spk_offset'])/1000.)
        else:
            raise KeyError('Spoke pattern {0:s} is not defined.'
                           .format(p['spk_pattern']))

        return w
    else:
        raise KeyError('Wheel template {0:s} is not defined.'.format(name))


def create_wheel_from_lm(l, mu):
    """Create a wheel with arbitrary values of lambda and mu.

    Since the reverse problem of generating a wheel with a specified lambda
    and mu is not unique, some constraints must be applied:

    - wheel dimensions are R=300m, d_hub=50mm, w_hub=25mm
    - radial spokes, E=200 GPa
    - number and diameter of spokes are chosen to ensure that Pn_lat = 2
    """

    wheel = BicycleWheel()

    wheel.hub = Hub(diam1=0.05, width1=0.025)

    wheel.rim = Rim(radius=0.3,
                    area=100e-6,
                    I11=1000e-12*(69./26.)*mu,
                    I22=1000e-12,
                    I33=1000e-12,
                    Iw=0.0e-12,
                    young_mod=69.0e9,
                    shear_mod=26.0e9)

    alpha = np.arctan(wheel.hub.width1 /
                      (wheel.rim.radius - wheel.hub.diam1/2))

    ls = np.sqrt(wheel.hub.width1**2 +
                 (wheel.rim.radius - wheel.hub.diam1/2)**2)

    # Calculate n_s by requiring that Pn_lat = 2.0 and round up to next even
    n_s = 1.5 * (8*np.pi/np.power(4, 0.25)) * np.power(l * (mu+1)/mu, 0.25)
    n_s = 2*int(np.ceil(n_s / 2))

    # Calculate spoke diameter based on requirement on lambda
    k_uu = l * wheel.rim.young_mod*wheel.rim.I22/(wheel.rim.radius**4)
    A_s = 2*np.pi*wheel.rim.radius*ls*k_uu / (n_s*200e9*np.sin(alpha)**2)
    d_s = np.sqrt(4*A_s/np.pi)

    # Create spokes
    wheel.lace_radial(n_spokes=n_s, diameter=d_s, young_mod=200e9)

    return wheel


def convert_to_paired(wheel):
    'Create a wheel with symmetric paired spokes with the same k_uu.'

    n_s = len(wheel.spokes)
    d_s = wheel.spokes[0].diameter / np.sqrt(2)
    E_s = wheel.spokes[0].EA / (np.pi*d_s**2/2)

    wheel.spokes = []

    theta = np.linspace(0.0, 2*np.pi, n_s+1)
    for t in theta[0:-1]:
        rim_pt = (wheel.rim.radius, t, 0.0)
        hub_pt_1 = (wheel.hub.radius1, t, wheel.hub.width1)
        hub_pt_2 = (wheel.hub.radius1, t, -wheel.hub.width2)

        s1 = Spoke(rim_pt, hub_pt_1, d_s, E_s)
        s2 = Spoke(rim_pt, hub_pt_2, d_s, E_s)

        wheel.spokes.append(s1)
        wheel.spokes.append(s2)

    return wheel
