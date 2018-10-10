#!/usr/bin/env python

'Post-process ABAQUS results file to a human-readable ASCII format'

import odbAccess
from abaqusConstants import *
import numpy as np
import sys
import os


def n_buckled(odb, frame_num, step_num=1):
    'Find number of buckled spokes.'

    part = odb.rootAssembly.instances['PART-1-1']
    regionSpokesAtRim = part.elementSets['ELSETSPOKESREF']

    step = odb.steps[odb.steps.keys()[step_num-1]]
    sf = step.frames[frame_num].fieldOutputs['SF']
    sf_spokes = sf.getSubset(region=regionSpokesAtRim)

    n_buck = 0

    for T in sf_spokes.values:
        if T.data[0] < 0.0:
            n_buck += 1

    return n_buck


odb_path = sys.argv[1]
odb_name = os.path.splitext(os.path.basename(odb_path))[0]

fname_Pd = odb_name + '_Pd.csv'
fname_shape = odb_name + '_shape.csv'

print odb_name

# Check if files already exist
overwrite = False
# Print load-deflection table
if not os.path.exists(fname_Pd) or overwrite:

    # Open ABAQUS output database
    odb = odbAccess.openOdb(odb_path)
    part = odb.rootAssembly.instances['PART-1-1']
    step = odb.steps[odb.steps.keys()[0]]

    with open(fname_Pd, 'w') as f_out:

        f_out.write('{0:s},{1:s},{2:s},{3:s},{4:s}\n'
                    .format('Time [s]', 'U2 [m]', 'RF2 [N]', 'U3 [m]',
                            'n_buckled'))

        for f_num in range(len(step.frames)):

            frame = step.frames[f_num]
            f_time = frame.frameValue

            # Displacements and reaction force
            U = frame.fieldOutputs['U']
            RF = frame.fieldOutputs['RF']
            U_load = U.getSubset(region=part.nodeSets['NSETLOAD'],
                                 position=NODAL).values[0]
            RF_load = RF.getSubset(region=part.nodeSets['NSETLOAD'],
                                   position=NODAL).values[0]

            n_buck = n_buckled(odb, f_num)

            f_out.write('{0:e},{1:e},{2:e},{3:e},{4:d}\n'
                        .format(f_time, float(U_load.data[1]),
                                float(RF_load.data[1]),
                                float(U_load.data[2]),
                                n_buck))
else:
    print('load-deflection data already exists')

# Print shapes file
N_shape_step = 10

if not os.path.exists(fname_shape) or overwrite:

    # Open ABAQUS output database
    odb = odbAccess.openOdb(odb_path)
    part = odb.rootAssembly.instances['PART-1-1']
    step = odb.steps[odb.steps.keys()[0]]

    with open(fname_shape, 'w') as f_out:

        reg = part.nodeSets['NSETSPOKENIP']

        # Print node labels
        f_out.write('    \tNode label->\t')
        f_out.write('\t'.join(['{0:12d}'.format(n.label) for n in reg.nodes])+'\n')

        # Print angular positions
        f_out.write('    \tAngle ->    ')
        for n in reg.nodes:
            x, y = (n.coordinates[0], n.coordinates[1])
            theta = np.arctan2(x + 1e-9, -y)
            if theta < 0:
                theta = 2*np.pi + theta
            f_out.write('\t{0: .5e}'.format(theta))
        f_out.write('\n')

        f_out.write('{0:4s}\t{0:11s}\n'.format('DOF', 'Time [s]'))

        def write_row(f_num, fieldVar, comp):
            field = step.frames[f_num].fieldOutputs[fieldVar]

            field_reg = field.getSubset(region=reg, position=NODAL)

            f_out.write('{0:2s} {1:d}\t'.format(fieldVar, comp))
            f_out.write('{0:.6e}\t'.format(step.frames[f_num].frameValue))
            f_out.write('\t'.join(['{0: .5e}'.format(float(val.data[comp]))
                                   for val in field_reg.values])+'\n')

        for f_num in range(0, len(step.frames), N_shape_step):
            write_row(f_num, 'U', 0)
            write_row(f_num, 'U', 1)
            write_row(f_num, 'U', 2)
            write_row(f_num, 'UR', 0)
            write_row(f_num, 'UR', 1)
            write_row(f_num, 'UR', 2)
else:
    print('shape data already exists')
