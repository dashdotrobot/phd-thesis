#!/usr/bin/env python

'Post-process ABAQUS results file to a human-readable ASCII format'

import odbAccess
from abaqusConstants import *
import numpy as np
import sys
import os


odb_path = sys.argv[1]
odb_name = os.path.splitext(os.path.basename(odb_path))[0]

fname_shape = odb_name + '_shape.csv'

print odb_name

# Check if files already exist
overwrite = False

# Print shapes file

if not os.path.exists(fname_shape) or overwrite:

    # Open ABAQUS output database
    odb = odbAccess.openOdb(odb_path)
    part = odb.rootAssembly.instances['PART-1-1']
    step = odb.steps[odb.steps.keys()[1]]

    with open(fname_shape, 'w') as f_out:

        reg = part.nodeSets['NSETSPOKENIP']

        # Print angular positions
        f_out.write('Angle,')
        for n in reg.nodes:
            x, y = (n.coordinates[0], n.coordinates[1])
            theta = np.arctan2(x + 1e-9, -y)
            if theta < 0:
                theta = 2*np.pi + theta
            f_out.write(',{0: .5e}'.format(theta))
        f_out.write('\n')

        f_out.write('DOF,Time [s],')
        f_out.write(','.join(['{0:12d}'.format(n.label) for n in reg.nodes])+'\n')

        def write_row(f_num, fieldVar, comp):
            field = step.frames[f_num].fieldOutputs[fieldVar]

            field_reg = field.getSubset(region=reg, position=NODAL)

            f_out.write('{0:s}{1:d},'.format(fieldVar, comp+1))
            f_out.write('{0:.6e},'.format(step.frames[f_num].frameValue))
            f_out.write(','.join(['{0: .5e}'.format(float(val.data[comp]))
                                  for val in field_reg.values])+'\n')

        # Write deformed shape
        write_row(1, 'U', 0)
        write_row(1, 'U', 1)
        write_row(1, 'U', 2)
        write_row(1, 'UR', 0)
        write_row(1, 'UR', 1)
        write_row(1, 'UR', 2)
else:
    print('shape data already exists')
