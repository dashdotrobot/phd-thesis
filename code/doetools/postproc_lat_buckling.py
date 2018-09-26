#!/usr/bin/env python

'Post-process ABAQUS results file to a human-readable ASCII format'

import odbAccess
from abaqusConstants import *
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

        f_out.write('{0:s},{1:s},{2:s},{3:s}\n'
                    .format('Time [s]', 'U3 [m]', 'RF3 [N]',
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

            f_out.write('{0:e},{1:e},{2:e},{3:d}\n'
                        .format(f_time, float(U_load.data[2]),
                                float(RF_load.data[2]),
                                n_buck))
else:
    print('load-deflection data already exists')
