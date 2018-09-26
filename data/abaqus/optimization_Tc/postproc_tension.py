#!/usr/bin/env python

'Post-process ABAQUS tension-buckling results file to a human-readable ASCII format'

import odbAccess
from abaqusConstants import *
import sys
import os


def get_node_field(odb, f_id, field_name, comp, n_id=0, nset=None):
    'Return the value of a field at a particular node.'

    part = odb.rootAssembly.instances['PART-1-1']
    step = odb.steps[odb.steps.keys()[0]]

    frame = step.frames[f_id]
    f_out = frame.fieldOutputs[field_name]

    if nset is not None:
        region = part.nodeSets[nset]
    else:
        region = part.nodes[n_id]

    d = f_out.getSubset(region=region,
                        position=NODAL).values[0].data[comp]
    return d


def N_avg(odb, f_id):
    'Return average rim compression in frame f_id.'

    # Get part information
    part = odb.rootAssembly.instances['PART-1-1']
    regionRim = part.elementSets['ELSETRIM']

    # Get step and field output
    step = odb.steps[odb.steps.keys()[0]]

    # Print table of spoke tension over time
    f = step.frames[f_id]
    sf = f.fieldOutputs['SF']

    sf_at_rim = sf.getSubset(region=regionRim,
                             position=INTEGRATION_POINT,
                             elementType='B31')

    N_avg = 0.0
    n = 0
    for N in sf_at_rim.values:
        N_avg = N_avg + N.data[0]
        n += 1
    N_avg = N_avg / n

    return N_avg


def T_avg(odb, f_id):
    'Returns the average spoke tension in frame f_id.'

    # Get part information
    part = odb.rootAssembly.instances['PART-1-1']
    regionSpokes = part.elementSets['ELSETSPOKES']

    # Get step and field output
    step = odb.steps[odb.steps.keys()[0]]

    # Print table of spoke tension over time
    f = step.frames[f_id]
    sf = f.fieldOutputs['SF']

    sf_at_rim = sf.getSubset(region=regionSpokes,
                             position=INTEGRATION_POINT,
                             elementType='B31')

    T_avg = 0.0
    n = 0
    for T in sf_at_rim.values:
        T_avg = T_avg + T.data[0]
        n += 1
    T_avg = T_avg / n

    return T_avg


def get_frame_temp(odb, f_id, elset='ELSETSPOKES'):
    'Get spoke temperature (1.0 = estimated buckling tension).'

    part = odb.rootAssembly.instances['PART-1-1']
    step = odb.steps[odb.steps.keys()[0]]
    frame = step.frames[f_id]

    if 'TEMP' in frame.fieldOutputs.keys():
        f_out = frame.fieldOutputs['TEMP']
        temp = f_out.getSubset(region=part.elementSets[elset],
                               position=INTEGRATION_POINT).values[0].data[0]
    else:
        temp = frame.frameValue

    return temp


verbose = False

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
        num_f = len(odb.steps[odb.steps.keys()[0]].frames)

        t = [get_frame_temp(odb, f) for f in range(num_f)]
        u = [get_node_field(odb, f, 'U', 2, n_id=0)
             for f in range(num_f)]
        T = [T_avg(odb, f) for f in range(num_f)]
        N = [N_avg(odb, f) for f in range(num_f)]

        f_out.write('{0:s},{1:s},{2:s},{3:s}\n'
                    .format('time', 'U3 [m]', 'Tension [N]',
                            'Compression [N]'))
        for tt, uu, TT, NN in zip(t, u, T, N):
            f_out.write('{0:e},{1:e},{2:e},{3:e}\n'.format(tt, uu, TT, NN))

        if verbose:
            print '{0:s},{1:s},{2:s},{3:s}\n'\
                .format('time', 'U3 [m]', 'Tension [N]', 'Compression [N]')
            for tt, uu, TT, NN in zip(t, u, T, NN):
                print '{0:e},{1:e},{2:e},{3:s}'.format(tt, uu, TT, NN)

else:
    print 'File {:s} already exists.'.format(fname_Pd)
