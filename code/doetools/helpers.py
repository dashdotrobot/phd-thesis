#!/usr/env/bin python

'Helper functions for ABAQUS doetools'

import numpy as np


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

        s1 = wheel.Spoke(rim_pt, hub_pt_1, d_s, E_s)
        s2 = wheel.Spoke(rim_pt, hub_pt_2, d_s, E_s)

        wheel.spokes.append(s1)
        wheel.spokes.append(s2)

    return wheel


def wheel_to_dict(wheel):
    'Create a dictionary of wheel parameters.'

    d = {'rim_radius': wheel.rim.radius,
         'rim_area': wheel.rim.area,
         'rim_I11': wheel.rim.I11,
         'rim_I22': wheel.rim.I22,
         'rim_I33': wheel.rim.I33,
         'rim_Iw': wheel.rim.Iw,
         'rim_young_mod': wheel.rim.young_mod,
         'rim_shear_mod': wheel.rim.shear_mod,
         'hub_width1': wheel.hub.width1,
         'hub_width2': wheel.hub.width2,
         'hub_diam1': wheel.hub.diam1,
         'hub_diam2': wheel.hub.diam2,
         'spk_num': len(wheel.spokes),
         'spk_diameter': wheel.spokes[0].diameter,
         'spk_young_mod': wheel.spokes[0].young_mod}

    return d
