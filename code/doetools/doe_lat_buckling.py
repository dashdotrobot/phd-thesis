#!/usr/bin/env python

'Design-of-experiment for lateral buckling numerical simulations.'

import pandas as pd
import numpy as np
from bikewheelcalc import ModeMatrix
from bikewheelcalc.continuum_analysis import calc_buckling_tension
from .abaqus_model import AbaqusModel
from .abaqus_doe import AbaqusDOE


abq_exp_preamble = """*Heading
** Non-linear post-buckling of pre-tensioned wheel under lateral force
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
**DAMPING, alpha={alpha:e}
*BOUNDARY
nsetHub, ENCASTRE\n"""

abq_exp_step = """**
*STEP, name=buckle
*DYNAMIC, explicit, fixed time incrementation
, {time:e}
*BOUNDARY, type=velocity
1, 3, , {vel:e}
1, 1, , 0.0
*OUTPUT, field, variable=PRESELECT, number interval=100
*ELEMENT OUTPUT, elset=elsetSpokes
SF
*ELEMENT OUTPUT, elset=elsetRim
SF
*OUTPUT, history, variable=PRESELECT
*END STEP\n"""

abq_riks_step = """** Riks post-buckling
*NSET, nset=nsetLoad
 1
*STEP, name=riks, inc=1000, nlgeom=YES
*STATIC, riks
 , , , 10.0, {sim_Pmax:e}, 1, 3, {sim_u3:e}
*BOUNDARY, fixed
nsetPreT, 1, 1
*CLOAD
 1, 3, 1.0
*OUTPUT, history, variable=PRESELECT
*OUTPUT, field, variable=PRESELECT
*ELEMENT OUTPUT, elset=elsetSpokes
SF
*ELEMENT OUTPUT, elset=elsetRim
SF
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
                '** Lateral collapse with {0:s} spokes\n'.format(opts['spk_eltype']) +
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


def write_abaqus_exp(fname, wheel, opts):
    'Write ABAQUS input file for collapse step.'

    if opts['spk_eltype'] == 'truss':
        am = AbaqusModel(wheel, n_spk=1)
    else:
        am = AbaqusModel(wheel, n_spk=10)

    # ABAQUS Explicit: dynamic collapse
    with open(fname, 'w') as f:
        f.write(abq_exp_preamble.format(alpha=opts['sim_alpha']))

        # Redefine rigid bodies
        f.write(am.write_rigid_ties())

        # Induce buckling
        f.write(abq_exp_step.format(time=opts['sim_time'],
                                    vel=opts['sim_u3']/opts['sim_time']))


def write_abaqus_riks(fname, wheel, opts):
    'Write (append) RIKS non-linear post-buckling step.'

    with open(fname, 'a') as f:
        f.write(abq_riks_step.format(sim_Pmax=opts['sim_Pmax'],
                                     sim_u3=opts['sim_u3']))


class LateralBucklingDOE(AbaqusDOE):
    'DOE for lateral buckling simulations'

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

            # Apply tension
            w.apply_tension(j['spk_T'])

            if j['sim_type'] == 'exp':

                # Write ABAQUS input files
                write_abaqus_pre(self.out_dir + '/' + jobname + '_preT.inp',
                                 w, j)
                write_abaqus_exp(self.out_dir + '/' + jobname + '_collapse.inp',
                                 w, j)

                # Write to batch file
                bf = batchfiles[b]
                bf.write('call abaqus interactive ')
                bf.write('job={0:s}_preT input={0:s}_preT.inp\n'
                         .format(jobname))

                bf.write('call abaqus interactive ')
                bf.write('job={0:s}_collapse input={0:s}_collapse.inp oldjob={0:s}_preT\n'
                         .format(jobname))

                # Next batchfile
                b = (b + 1) % N_batches

                # Write entry to postprocess script
                pp_script.write('call abaqus python postproc_lat_buckling.py ' +
                                '{0:s}_collapse.odb\n'.format(jobname))

            elif j['sim_type'] == 'riks':

                # Write ABAQUS input file
                write_abaqus_pre(self.out_dir + '/' + jobname + '_riks.inp',
                                 w, j)
                write_abaqus_riks(self.out_dir + '/' + jobname + '_riks.inp',
                                  w, j)

                # Write to batch file
                bf = batchfiles[b]
                bf.write('call abaqus interactive ')
                bf.write('job={0:s}_riks input={0:s}_riks.inp\n'
                         .format(jobname))

                # Next batchfile
                b = (b + 1) % N_batches

                # Write entry to postprocess script
                pp_script.write('call abaqus python postproc_lat_riks.py ' +
                                '{0:s}_riks.odb\n'.format(jobname))

        # Close batch files
        for f in batchfiles:
            f.close()

    def extract_results(self):
        'Extract results and calculated quantities from ABAQUS output'

        def spoke_buckle_load(pd_data):
            'Get load at which spokes start to buckle'

            return pd_data.RF3[np.nonzero(pd_data.n_buckled)[0][0]]

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
                if j['sim_type'] == 'exp':
                    pd_fname = self.out_dir + '/' + j.name + '_collapse_Pd.csv'
                elif j['sim_type'] == 'riks':
                    pd_fname = self.out_dir + '/' + j.name + '_riks_Pd.csv'

                pd_data = pd.read_csv(pd_fname)
                pd_data.columns = ['Time', 'U3', 'RF3', 'n_buckled']

                # Stiffness
                w.apply_tension(0.01)
                K_lat_0 = mm.calc_lat_stiff(smeared_spokes=True, coupling=False)
                K_rad_0 = mm.calc_rad_stiff(smeared_spokes=False, coupling=False)

                self.db.at[i, 'K_lat_0'] = K_lat_0
                self.db.at[i, 'K_rad_0'] = K_rad_0

                # Critical loads
                self.db.at[i, 'Pc_max'] = max(pd_data.RF3)
                self.db.at[i, 'Pc_spk'] = spoke_buckle_load(pd_data)

            except Exception as e:
                print('Error on {0:s}: {1:s}'.format(j.name, str(e)))
                continue

    def __init__(self, out_dir, db_file=None, opts={}):
        'Create design-of-experiment'

        opts_default = {'sim_type': 'exp',  # 'exp' or 'riks'
                        'sim_time': 1.0,    # [s] Simulation time
                        'sim_u3': 0.04,     # [w] Lateral displacement
                        'spk_T': 0.0,       # [N] Spoke tension
                        'sim_alpha': 1.e3   # Material damping
                        }

        # Call parent constructor
        opts_default.update(opts)
        AbaqusDOE.__init__(self, out_dir, db_file, opts_default)
