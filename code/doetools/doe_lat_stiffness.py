#!/usr/bin/env python

'Design-of-experiment for tension-buckling simulations'

import pandas as pd
import numpy as np
from numpy import pi, cos
from .abaqus_model import AbaqusModel
from .abaqus_doe import AbaqusDOE


class LateralStiffnessDOE(AbaqusDOE):

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

            # Apply pretension
            w.apply_tension(j['spk_T'])

            if j['spk_eltype'] == 'truss':
                am = AbaqusModel(w, n_spk=1)
            else:
                am = AbaqusModel(w, n_spk=10)

            # Write ABAQUS input file
            with open(self.out_dir + '/' + jobname + '.inp', 'w') as f:
                f.write('*Heading\n' +
                        '** Tension buckling with {0:s} spokes\n'
                        .format(j['spk_eltype']) +
                        '**\n')

                f.write(am.write_heading('Nodes'))
                f.write(am.write_rim_nodes())
                f.write(am.write_spoke_nodes())
                f.write(am.write_pretension_nodes())
                f.write(am.write_rigid_ties())

                f.write(am.write_heading('Elements'))
                f.write(am.write_rim_elems())

                if j['spk_eltype'] == 'truss':
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
                        '*BOUNDARY, op=new\n' +
                        'nsetHub, ENCASTRE\n' +
                        '*CLOAD\n')

                for i, s in enumerate(w.spokes):
                    f.write(' {:5d}, 1, {:e}\n'.format(99000 + i+1, s.tension))

                f.write('*OUTPUT, field, variable=PRESELECT\n' +
                        '*ELEMENT OUTPUT, elset=elsetSpokes\nSF\n' +
                        '*ELEMENT OUTPUT, elset=elsetRim\nSF\n' +
                        '*OUTPUT, history, variable=PRESELECT\n' +
                        '*END STEP\n')

                f.write(am.write_heading('Lateral Stiffness'))
                f.write('*STEP, name=K_lat, perturbation\n' +
                        '*STATIC\n' +
                        '*BOUNDARY, fixed\nnsetPreT, 1\n' +
                        '*CLOAD\n1, 3, 1\n' +
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
            pp_script.write('call abaqus python postproc_lat_stiffness.py ' +
                            '{0:s}.odb\n'.format(jobname))

        # Close batch files
        for f in batchfiles:
            f.close()

    def extract_results(self):
        'Extract results and calculated quantities from ABAQUS output'

        for i in self.db.index:
            print('.', end='')

            j = self.db.loc[i]

            try:
                # Get deflection at load point
                shape = pd.read_csv(self.out_dir + '/' + j.name +
                                    '_shape.csv',
                                    header=[1])

                self.db.at[i, 'K_lat_abq'] = 1.0 /\
                    shape[shape['DOF'] == 'U3'].iloc[0][2]

            except Exception as e:
                print('Error on {0:s}: {1:s}'.format(j.name, e))
                continue

    def __init__(self, out_dir, db_file=None, opts={}):
        'Create design-of-experiment for tension buckling'

        opts_default = {'spk_T': 0.}

        # Call parent constructor
        opts_default.update(opts)
        AbaqusDOE.__init__(self, out_dir, db_file, opts_default)
