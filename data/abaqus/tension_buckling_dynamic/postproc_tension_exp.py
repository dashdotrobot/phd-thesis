#!/usr/bin/env python

'Post-process ABAQUS tension-buckling results file to a human-readable ASCII format'

import odbAccess
from abaqusConstants import *
import sys
import os


def get_node_field(odb, step_id, f_id, field_name, comp, n_id=0, nset=None):
    'Return the value of a field at a particular node.'

    part = odb.rootAssembly.instances['PART-1-1']
    step = odb.steps[step_id]

    frame = step.frames[f_id]
    f_out = frame.fieldOutputs[field_name]

    if nset is not None:
        region = part.nodeSets[nset]
    else:
        region = part.nodes[n_id]

    d = f_out.getSubset(region=region,
                        position=NODAL).values[0].data[comp]
    return d


def get_N_avg(odb, step_id, f_id):
    'Return average rim compression in frame f_id.'

    # Get part information
    part = odb.rootAssembly.instances['PART-1-1']
    regionRim = part.elementSets['ELSETRIM']

    # Get step and field output
    step = odb.steps[step_id]

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


def get_T_avg(odb, step_id, f_id):
    'Returns the average spoke tension in frame f_id.'

    # Get part information
    part = odb.rootAssembly.instances['PART-1-1']
    regionSpokes = part.elementSets['ELSETSPOKESREF']

    # Get step and field output
    step = odb.steps[step_id]

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


def get_frame_temp(odb, step_id, f_id, elset='ELSETSPOKESREF'):
    'Get spoke temperature (1.0 = estimated buckling tension).'

    part = odb.rootAssembly.instances['PART-1-1']
    step = odb.steps[step_id]
    frame = step.frames[f_id]

    if 'TEMP' in frame.fieldOutputs.keys():
        f_out = frame.fieldOutputs['TEMP']
        temp = f_out.getSubset(region=part.elementSets[elset],
                               position=INTEGRATION_POINT).values[0].data
    else:
        temp = frame.frameValue

    return temp


verbose = False

odb_path = sys.argv[1]
odb_name = os.path.splitext(os.path.basename(odb_path))[0]

fname_Pd = odb_name + '_Pd.csv'
fname_SE = odb_name + '_SE.csv'

# Check if files already exist
overwrite = False

# Print load-deflection table
if not os.path.exists(fname_Pd) or overwrite:

    # Setup result vectors
    s = []
    time = []
    temp = []
    u3 = []
    Tavg = []
    Navg = []

    # Open ABAQUS output database
    odb = odbAccess.openOdb(odb_path)
    part = odb.rootAssembly.instances['PART-1-1']

    s_num = 0
    for s_id in odb.steps.keys():
        step = odb.steps[s_id]

        n_frames = len(step.frames)

        s.extend([s_num]*n_frames)

        time.extend([step.frames[f_id].frameValue for f_id in range(n_frames)])

        temp.extend([get_frame_temp(odb, s_id, f_id)
                     for f_id in range(n_frames)])
        u3.extend([get_node_field(odb, s_id, f_id, 'U', 2, n_id=0)
                   for f in range(n_frames)])
        Tavg.extend([get_T_avg(odb, s_id, f_id)
                     for f_id in range(n_frames)])
        Navg.extend([get_N_avg(odb, s_id, f_id)
                     for f_id in range(n_frames)])

        s_num += 1

    with open(fname_Pd, 'w') as f_out:

        f_out.write('{0:s},{1:s},{2:s},{3:s},{4:s},{5:s}\n'
                    .format('Step', 'Time [s]', 'Temperature', 'U3 [m]',
                            'Tension [N]', 'Compression [N]'))
        for ss, tm, tp, u, T, N in zip(s, time, temp, u3, Tavg, Navg):
            f_out.write('{0:d},{1:e},{2:e},{3:e},{4:e},{5:e}\n'
                        .format(ss, tm, tp, u, T, N))

        if verbose:
            print '{0:12s} {1:12s} {2:12s} {3:12s} {4:12s}\n'\
                .format('Time [s]', 'Temperature', 'U3 [m]', 'Tension [N]',
                        'Compression [N]')
            for tm, tp, u, T, N in zip(time, temp, u3, Tavg, Navg):
                print '{0:6e} {1:6e} {2:6e} {3:6e} {4:6e}'\
                    .format(tm, tp, u, T, N)

else:
    print 'File {:s} already exists.'.format(fname_Pd)


if not os.path.exists(fname_SE) or overwrite:

    s = []
    time = []
    strain_energy = []

    # Open ABAQUS output database
    odb = odbAccess.openOdb(odb_path)
    part = odb.rootAssembly.instances['PART-1-1']

    s_num = 0
    for s_id in odb.steps.keys():
        step = odb.steps[s_id]

        reg = step.historyRegions[step.historyRegions.keys()[0]]
        t, se = zip(*reg.historyOutputs['ALLSE'].data)

        s.extend([s_num]*len(t))
        time.extend(t)
        strain_energy.extend(se)

        s_num += 1

    with open(fname_SE, 'w') as f_SE:
        f_SE.write('Step,Time [s],Strain Energy [J]\n')

        for st, t, se in zip(s, time, strain_energy):
            f_SE.write('{0:d},{1:e},{2:e}\n'.format(st, t, se))

else:
    print 'File {:s} already exists.'.format(fname_SE)
