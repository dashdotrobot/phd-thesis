#!/usr/bin/env python

'Base class for ABAQUS numerical experiments.'

import pandas as pd
import pickle
from bikewheelcalc import BicycleWheel, Hub, Rim
from .helpers import convert_to_paired, wheel_to_dict
import hashlib


class AbaqusDOE:
    'Base class for ABAQUS numerical experiments'

    def wheel_from_row(self, j):
        'Create a BicycleWheel object from a job row'

        # Try to un-pickle wheel object
        if 'wheel_pickle' in j:
            with open(self.out_dir + '/' + j['wheel_pickle'], 'rb') as p:
                w = pickle.load(p)

        else:
            w = BicycleWheel()
            w.hub = Hub(diam1=j['hub_diam1'], width1=j['hub_width1'])
            w.rim = Rim(radius=j['rim_radius'],
                        area=j['rim_area'],
                        I11=j['rim_I11'],
                        I22=j['rim_I22'],
                        I33=j['rim_I33'],
                        Iw=j['rim_Iw'],
                        young_mod=j['rim_young_mod'],
                        shear_mod=j['rim_shear_mod'])

            w.lace_radial(n_spokes=int(j['spk_num']),
                          diameter=j['spk_diameter'],
                          young_mod=j['spk_young_mod'])

            if 'spk_paired' in j.keys() and j['spk_paired']:
                w = convert_to_paired(w)

        return w

    def add_experiment(self, wheel, opts={}):
        'Add experiment to the current design-of-experiment'

        # Set job options
        job_opts = self.opts_global.copy()
        job_opts.update(opts)

        # Add to DOE database
        job_opts.update(wheel_to_dict(wheel))

        job_series = pd.Series(job_opts)

        if 'jobname' in job_opts:
            job_series.name = job_opts['jobname']
        else:
            # first 6 digits of hash
            job_series.name = hashlib.md5(str(job_series)).hexdigest()[0:6]

        # Pickle wheel object
        job_series['wheel_pickle'] = '{0:s}.pkl'.format(job_series.name)
        with open(self.out_dir + '/' + job_series['wheel_pickle'], 'wb') as p:
            pickle.dump(wheel, p)

        self.db = self.db.append(job_series)

    def to_csv(self):
        'Write DOE experiment database to CSV'

        self.db.to_csv(self.out_dir + '/_doe_db.csv')

    def read_csv(self, db_file):
        'Load DOE data from CSV'

        def convert_rim_perturb(x):
            if isinstance(x, str):
                return map(float, x.strip('[]').split(','))
            else:
                return None

        self.db = pd.read_csv(db_file, index_col=0)

        # Convert rim_perturb to list or None
        self.db['rim_perturb'] = \
            self.db['rim_perturb'].map(convert_rim_perturb)

    def __init__(self, out_dir, db_file=None, opts={}):
        'Create design-of-experiment'

        # Output directory
        self.out_dir = out_dir

        # Default job options
        self.opts_global = {'rim_perturb': None,   # Perturbation coeffs
                            'spk_eltype': 'beam',  # 'truss' or 'beam'
                            }

        self.opts_global.update(opts)

        if db_file is not None:  # Read existing database from CSV
            self.read_csv(db_file)
            self.opts_global = {}
        else:
            self.db = pd.DataFrame()
