#!/usr/bin/env python

'Design-of-experiment for tension-buckling simulations'

import pandas as pd
import numpy as np
from numpy import pi, cos
from bikewheelcalc.continuum_analysis import calc_buckling_tension
from .abaqus_model import AbaqusModel
from .abaqus_doe import AbaqusDOE


class TensionBucklingDOE(AbaqusDOE):

    def write_input_files(self, N_batches=1):
        # Open list of batch files
        batchfiles = []
        for b in range(N_batches):
            fname = '{0:s}/_run_{1:d}.bat'.format(self.out_dir, b)
            batchfiles.append(open(fname, 'w'))

        pp_script = open('{0:s}/_postproc.bat'.format(self.out_dir), 'w')

        b = 0
        for jobname, j in self.db.iterrows():
            w = self.wheel_from_row(j)

            R = w.rim.radius
            ns = len(w.spokes)

            # Estimate buckling tension
            Tc, nc = calc_buckling_tension(w)
            vc = R*ns*Tc/(2*pi)/(w.rim.young_mod*w.rim.area)

            w.apply_tension(Tc)

            # Estimate temperature to cause buckling
            deltaT = -j['sim_temp_factor']

            alpha1 = (w.spokes[0].tension/w.spokes[0].EA +
                      vc*w.spokes[0].n[1]/w.spokes[0].length)
            alpha2 = (w.spokes[1].tension/w.spokes[1].EA +
                      vc*w.spokes[1].n[1]/w.spokes[1].length)

            if j['rim_perturb'] is None:
                perturb_a = [0.00e-3, 0.01e-3, 0.01e-3]
            else:
                perturb_a = j['rim_perturb']

            if j['spk_eltype'] == 'truss':
                am = AbaqusModel(w, n_spk=1)
            else:
                am = AbaqusModel(w, n_spk=10)

            def perturb_rim(theta):
                u = 0.0
                for n in range(len(perturb_a)):
                    u = u + perturb_a[n]*cos((n+1)*theta)
                return u

            # Write ABAQUS input file
            with open(self.out_dir + '/' + jobname + '.inp', 'w') as f:
                f.write('*Heading\n' +
                        '** Tension buckling with {0:s} spokes\n'.format(j['spk_eltype']) +
                        '**\n')

                f.write(am.write_heading('Nodes'))
                f.write(am.write_rim_nodes(f_perturb=perturb_rim))
                f.write(am.write_spoke_nodes())
                f.write(am.write_rigid_ties())

                f.write(am.write_heading('Elements'))
                f.write(am.write_rim_elems())

                if j['spk_eltype'] == 'truss':
                    f.write(am.write_spoke_elems(eltype='T3D2'))
                else:
                    f.write(am.write_spoke_elems())

                f.write(am.write_heading('Sections'))
                f.write(am.write_beam_sections(alpha1=alpha1, alpha2=alpha2))

                f.write(am.write_heading('Boundary conditions'))
                f.write(am.write_bc_fix_hub())

                f.write(am.write_heading('Apply tension'))
                f.write('*STEP, name=buckle, nlgeom=yes\n' +
                        '*STATIC\n' +
                        ' 1.0e-3, , , 2.0e-2\n' +
                        '*TEMPERATURE\n' +
                        'nsetSpokes, {temp:e}\n'.format(temp=deltaT) +
                        '*OUTPUT, field, variable=PRESELECT\n' +
                        '*ELEMENT OUTPUT, elset=elsetSpokes\nSF\n' +
                        '*ELEMENT OUTPUT, elset=elsetRim\nSF\n' +
                        '*OUTPUT, history, variable=PRESELECT\n' +
                        '*END STEP\n')

            # Write to batch file
            bf = batchfiles[b]
            bf.write('call abaqus interactive ')
            bf.write('job={0:s} input={0:s}.inp\n'.format(jobname))

            # Next batchfile
            b = (b + 1) % N_batches

            # Write entry to postprocess script
            pp_script.write('call abaqus python postproc_tension.py ' +
                            '{0:s}.odb\n'.format(jobname))

        # Close batch files
        for f in batchfiles:
            f.close()

    def extract_results(self):
        'Extract results and calculated quantities from ABAQUS output'

        def Tc_Southwell(pd_data, T_max):
            'Calculate Tc using Southwells method'

            u = pd_data['U3 [m]'].values[1:]
            T = pd_data['Tension [N]'].values[1:]

            N = int(0.8*np.where(T >= T_max)[0][0])

            p = np.polyfit(u[:N], u[:N]/T[:N], 2)

            return 1. / (2*p[0]*u[N] + p[1])

        def Tc_nonlin(pd_data, tol=0.02):
            'Calculate Tc from departure from linearity'

            # Fit straight line to first 10% of data
            x, y = (pd_data['time'][:-1], pd_data['Tension [N]'][:-1])
            N = int(np.max([10, np.ceil(0.10 * len(x))]))
            pf = np.polyfit(x[:N], y[:N], 2)

            err = (np.polyval(pf, x) - y) / np.max(y)
            return y[np.argmax(np.abs(err) > tol)]

        for i in self.db.index:
            print('.', end='')

            j = self.db.loc[i]
            w = self.wheel_from_row(j)

            try:
                # Buckling tension
                Tc, nc = calc_buckling_tension(w)
                self.db.at[i, 'Tc'] = Tc
                self.db.at[i, 'nc'] = nc

                pd_data = pd.read_csv(self.out_dir + '/' + j.name + '_Pd.csv')
                self.db.at[i, 'Tc_nonlin'] = Tc_nonlin(pd_data)
                self.db.at[i, 'Tc_south'] = \
                    Tc_Southwell(pd_data, self.db.at[i, 'Tc_nonlin'])

            except Exception as e:
                print('Error on {0:s}: {1:s}'.format(j.name, str(e)))
                continue

    def __init__(self, out_dir, db_file=None, opts={}):
        'Create design-of-experiment for tension buckling'

        opts_default = {'sim_temp_factor': 1.5}

        # Call parent constructor
        opts_default.update(opts)
        AbaqusDOE.__init__(self, out_dir, db_file, opts_default)


if False:
    import wheel_library
    doe = TensionBucklingDOE(out_dir='res_test_tb')

    wheel = wheel_library.create_predefined_wheel('vintage_rd')
    doe.add_experiment(wheel)
    doe.add_experiment(wheel, opts={'rim_perturb': [0, 1, 2]})
    doe.add_experiment(wheel, opts={'sim_temp_factor': 2.0})
    doe.add_experiment(wheel, opts={'sim_time': 1.0, 'spk_eltype': 'truss'})

    doe.write_input_files(N_batches=2)
    doe.to_csv()

    doe2 = TensionBucklingDOE(out_dir='res_test_tb2',
                              db_file='res_test_tb/_doe_db.csv')

    doe2.write_input_files()
    doe2.to_csv()
