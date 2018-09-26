#!/usr/bin/env python

'Design-of-experiment for radial buckling numerical simulations.'

import pandas as pd
import numpy as np
from bikewheelcalc import ModeMatrix
from bikewheelcalc.continuum_analysis import calc_buckling_tension
from .abaqus_model import AbaqusModel
from .abaqus_doe import AbaqusDOE


abq_exp_preamble = """*Heading
** Non-linear post-buckling of pre-tensioned wheel under radial force
** EXPLICIT: dynamic post-buckling analysis
*Preprint, echo=NO, model=NO, history=NO, contact=NO
*IMPORT, update=yes
elsetRim, elsetSpokes,
*IMPORT elset
elsetRim, elsetSpokes, elsetSpokes1, elsetSpokes2, elsetSpokesRef
*IMPORT nset
nsetRim, nsetSpokes, nsetSpokeNip, nsetHub, nsetPreT
*NSET, nset=nsetLoad
 1
**
**MATERIAL, name=steel
**ELASTIC
** 2.1e+11, 0.33
**DENSITY
** 8050.0
**DAMPING, alpha=1000.0
*BOUNDARY
nsetHub, ENCASTRE\n"""

abq_exp_step = """**
*STEP, name=buckle
*DYNAMIC, explicit, fixed time incrementation
, {time:e}
*BOUNDARY, type=velocity
1, 2, , {vel:e}
1, 1, , 0.0
*OUTPUT, field, variable=PRESELECT, number interval=100
*ELEMENT OUTPUT, elset=elsetSpokes
SF
*ELEMENT OUTPUT, elset=elsetRim
SF
*OUTPUT, history, variable=PRESELECT
*END STEP\n"""


def write_abaqus_pre(fname, wheel, opts):
    'Write ABAQUS input file for pretension step.'

    if 'perturb_a' not in opts:
        perturb_a = [0.01e-3, 0.01e-3, 0.01e-3]
    else:
        perturb_a = opts['perturb_a']

    def perturb_rim(theta):
        u = 0.0
        for n in range(len(perturb_a)):
            u = u + perturb_a[n]*np.cos((n+1)*theta)
        return u

    if opts['spk_eltype'] == 'truss':
        am = AbaqusModel(wheel, n_spk=1)
    else:
        am = AbaqusModel(wheel, n_spk=10)

    # ABAQUS Standard: apply spoke tension
    with open(fname, 'w') as f:

        f.write('*Heading\n' +
                '** Radial collapse with {0:s} spokes\n'.format(opts['spk_eltype']) +
                '**\n')

        f.write(am.write_heading('Nodes'))
        f.write(am.write_rim_nodes(f_perturb=perturb_rim))
        f.write(am.write_spoke_nodes())
        f.write(am.write_pretension_nodes())
        f.write(am.write_rigid_ties())

        f.write(am.write_heading('Elements'))
        f.write(am.write_rim_elems())

        if opts['spk_eltype'] == 'truss':
            f.write(am.write_spoke_elems(eltype='T3D2'))
        else:
            f.write(am.write_spoke_elems())

        f.write(am.write_heading('Sections'))
        f.write(am.write_pretension_section())
        f.write(am.write_beam_sections())

        f.write(am.write_heading('Boundary conditions'))
        f.write(am.write_bc_fix_hub())

        # Apply pretension
        f.write(am.write_heading('Pre-load spokes'))

        f.write('*STEP, name=preT, nlgeom=YES\n' +
                '*STATIC\n' +
                '1., 1., 1e-05, 1.\n' +
                '*CLOAD\n')

        for i, s in enumerate(wheel.spokes):
            f.write(' {:5d}, 1, {:e}\n'.format(99000 + i+1, s.tension))

        f.write('*OUTPUT, field, variable=PRESELECT\n' +
                '*ELEMENT OUTPUT, elset=elsetSpokes\nSF\n' +
                '*ELEMENT OUTPUT, elset=elsetRim\nSF\n' +
                '*RESTART, write, number interval=1\n'
                '*END STEP\n')


def write_abaqus_exp(fname, wheel, opts, u_rad):
    'Write ABAQUS input file for collapse step.'

    if opts['spk_eltype'] == 'truss':
        am = AbaqusModel(wheel, n_spk=1)
    else:
        am = AbaqusModel(wheel, n_spk=10)

    # ABAQUS Explicit: dynamic collapse
    with open(fname, 'w') as f:
        f.write(abq_exp_preamble)

        # Redefine rigid bodies
        f.write(am.write_rigid_ties())

        # Induce buckling
        f.write(abq_exp_step.format(time=opts['sim_time'],
                                    vel=u_rad/opts['sim_time']))


class RadialBucklingDOE(AbaqusDOE):
    'DOE for radial buckling simulations'

    def write_input_files(self, N_batches=1):
        'Write ABAQUS input files and batch script(s)'

        # Open list of batch files
        batchfiles = []
        for b in range(N_batches):
            fname = '{0:s}/_run_{1:d}.bat'.format(self.out_dir, b)
            batchfiles.append(open(fname, 'w'))

        pp_script = open('{0:s}/_postproc.bat'.format(self.out_dir), 'w')

        b = 0
        for jobname, j in self.db.iterrows():
            w = self.wheel_from_row(j)
            mm = ModeMatrix(w)

            # Apply tension
            if j['spk_T'] == 0.:
                w.apply_tension(0.001)
            else:
                w.apply_tension(j['spk_T'])

            # Choose total radial displacement
            if j['sim_u2'] == 'auto':
                # Estimate radial displacement to failure for truss wheel
                Pc_est = mm.calc_lat_stiff() * w.rim.radius
                u_rad = 1.5*Pc_est / mm.calc_rad_stiff()
            else:
                u_rad = j['sim_u2']

            # Write ABAQUS input files
            write_abaqus_pre(self.out_dir + '/' + jobname + '_preT.inp',
                             w, j)
            write_abaqus_exp(self.out_dir + '/' + jobname + '_collapse.inp',
                             w, j, u_rad)

            # Write to batch file
            bf = batchfiles[b]
            bf.write('call abaqus interactive ')
            bf.write('job={0:s}_preT input={0:s}_preT.inp\n'.format(jobname))

            bf.write('call abaqus interactive ')
            bf.write('job={0:s}_collapse input={0:s}_collapse.inp oldjob={0:s}_preT\n'.format(jobname))

            # Next batchfile
            b = (b + 1) % N_batches

            # Write entry to postprocess script
            pp_script.write('call abaqus python postproc_rad_buckling.py ' +
                            '{0:s}_collapse.odb\n'.format(jobname))

        # Close batch files
        for f in batchfiles:
            f.close()

    def extract_results(self):
        'Extract results and calculated quantities from ABAQUS output'

        def buckling_load_southwell(pd_data):
            'Get critical buckling load from Southwell plot'
            pass

        def buckling_load_nonlin(pd_data):
            'Get critical buckling from departure from linearity'

            # Fit straight line to first 5% of data
            N = int(np.ceil(0.05*len(pd_data)))
            pf = np.polyfit(pd_data['U2'][:N], pd_data['RF2'][:N], 1)

            err = (np.polyval(pf, pd_data['U2']) - pd_data['RF2']) /\
                np.mean(pd_data['RF2'])
            return pd_data['RF2'][np.argmax(np.abs(err) > 0.02)]

        def calc_T_nb(w, K_lat_0, K_rad_0, Tc, form=1):
            'Calculate minimum spoke tension for no spoke buckling'

            R = w.rim.radius
            EA = w.spokes[0].EA
            alpha = w.spokes[0].alpha
            ls = w.spokes[0].length

            if form == 2:
                # Quadratic approximation for K_lat(T)
                x = K_rad_0*Tc*ls/(2*EA*K_lat_0*R*np.cos(alpha))
                return np.sqrt(1+x**2) - x
            else:
                # Linear approximation for K_lat(T)
                return 1.0/(1 + K_rad_0*ls*Tc/(K_lat_0*R*EA*np.cos(alpha)))

        def calc_P_nb(w, K_lat_0, K_rad_0, Tc, form=1):
            'Critical load at the critical no-buckling tension'

            R = w.rim.radius
            EA = w.spokes[0].EA
            alpha = w.spokes[0].alpha
            ls = w.spokes[0].length

            if form == 2:
                # Quadratic approximation for K_lat(T)
                x = K_rad_0*Tc*ls/(2*EA*K_lat_0*R*np.cos(alpha))
                T_nb = np.sqrt(1+x**2) - x
                return K_lat_0*R*(1-T_nb**2)
            else:
                # Linear approximation for K_lat(T)
                return K_lat_0*R / (1 + EA/Tc*K_lat_0/K_rad_0)

        for i in self.db.index:
            print('.', end='')

            j = self.db.loc[i]
            w = self.wheel_from_row(j)
            mm = ModeMatrix(w)

            try:
                # Buckling tension
                Tc, nc = calc_buckling_tension(w)
                self.db.at[i, 'Tc'] = Tc
                self.db.at[i, 'nc'] = nc

                # Get load-deflection data
                pd_data = pd.read_csv(self.out_dir + '/' + j.name +
                                      '_collapse_Pd.csv')
                pd_data.columns = ['Time', 'U2', 'RF2', 'U3', 'n_buckled']

                # Stiffness at zero tension
                w.apply_tension(0.01)
                K_lat_0 = mm.calc_lat_stiff(smeared_spokes=True, coupling=False)
                K_rad_0 = mm.calc_rad_stiff(smeared_spokes=False, coupling=False)

                w.apply_tension(j['spk_T'])
                K_lat = mm.calc_lat_stiff(smeared_spokes=True, coupling=False)

                self.db.at[i, 'K_lat_0'] = K_lat_0
                self.db.at[i, 'K_rad_0'] = K_rad_0
                self.db.at[i, 'K_lat'] = K_lat

                # Radial buckling load
                self.db.at[i, 'Pc_max'] = max(pd_data.RF2)
                self.db.at[i, 'Pc_nonlin'] = buckling_load_nonlin(pd_data)

                self.db.at[i, 'T_nb'] = calc_T_nb(w, K_lat_0, K_rad_0, Tc)
                self.db.at[i, 'Pc_nb'] = calc_P_nb(w, K_lat_0, K_rad_0, Tc)
                self.db.at[i, 'T_nb_2'] = calc_T_nb(w, K_lat_0, K_rad_0, Tc,
                                                    form=2)
                self.db.at[i, 'Pc_nb_2'] = calc_P_nb(w, K_lat_0, K_rad_0, Tc,
                                                     form=2)

            except Exception as e:
                print('Error on {0:s}: {1:s}'.format(j.name, str(e)))
                continue

    def __init__(self, out_dir, db_file=None, opts={}):
        'Create design-of-experiment'

        opts_default = {'sim_time': 1.0,  # Simulation time
                        'sim_u2': 0.01,   # Radial displacement
                        'spk_T': 0.0      # Average spoke tension
                        }

        # Call parent constructor
        opts_default.update(opts)
        AbaqusDOE.__init__(self, out_dir, db_file, opts_default)
