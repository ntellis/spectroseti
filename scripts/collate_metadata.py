import spectroseti.apf as apf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import seaborn as sb
import spectroseti.apf as apf
import spectroseti.apfdefinitions as apfdefs
from os import listdir, mkdir
import spectroseti.utilities as util
import csv
import scipy.signal as sg
import pickle

import pandas as pd
from xarray import Dataset


run_name = 'bkjRunJun16'
folder_name = apfdefs.laser_search_run_dir + run_name + '/'

all_reduced = listdir(folder_name)
files_to_collate = [fn for fn in all_reduced if fn[-2:] == '.p']


d = dict()
all_devs = []
for f in files_to_collate:
    pickle_off = open(folder_name+f, "rb")
    emp = pickle.load(pickle_off)
    name = f[1:8]
    flat_dict = dict()
    flat_dict.update(emp['meta'])
    devs = emp['devs']
    meta = emp['meta']
    order_medians = emp['order_medians']
    color_index = order_medians[61]/order_medians[21]
    # devs_list = []
    run = meta['run']
    obs = meta['obs']
    for i, dev in enumerate(devs):
        counts_per_mad = dev['dev'][4]
        all_devs.append(
            {
                'run': run,
                'obs': obs,
                'cosmic_reject_value': dev['cosmic_reject_value'],
                'order': dev['dev'][0],
                'num_pixels': dev['dev'][1],
                'start_pixel': dev['dev'][2],
                'continuum_val': dev['dev'][3],
                'MAD': dev['dev'][4],
                'threshold': dev['dev'][3] + dev['dev'][4] *  meta['number_mads'],
                'peak_pixel_value': (dev['dev'][8] - dev['dev'][3])/counts_per_mad,
                'mean_deviant_pixel_value': (dev['dev'][5] - dev['dev'][3])/counts_per_mad,
                'median_deviant_pixel_value': (dev['dev'][6] - dev['dev'][3])/counts_per_mad,
                'central_wav': dev['dev'][7],
                'intensity': (dev['dev'][5] - dev['dev'][3]) * dev['dev'][1] / counts_per_mad,
                'color_index': color_index,
                'target_name': meta['target_name'],
                'exposure_time': meta['exposure_time'],
                'RA': meta['RA'],
                'DEC': meta['DEC'],
                'HA': meta['HA'],
                'AZ': meta['AZ'],
                'reduced_filename': 'r%(run)s.%(obs)s.fits' % locals(),
                'raw_filename': 'ucb-%(run)s%(obs)s.fits' % locals(),
                'img_prefix': 'r%(run)s.%(obs)s_dev%(i)s' % locals()
            }
        )
    # Now we should make a .csv of all of these data
    keys = ['run',
            'obs',
            'cosmic_reject_value',
            'order',
            'num_pixels',
            'start_pixel',
            'continuum_val',
            'MAD',
            'threshold',
            'peak_pixel_value',
            'mean_deviant_pixel_value',
            'median_deviant_pixel_value',
            'central_wav',
            'intensity',
            'color_index',
            'target_name',
            'exposure_time',
            'RA',
            'DEC',
            'HA',
            'AZ',
            'reduced_filename',
            'raw_filename',
            'img_prefix']
    filename = apfdefs.laser_search_run_dir + '/' + run_name +'_metadata.csv'
    with open(filename, 'wb') as output_file:
        writer = csv.DictWriter(
            output_file, fieldnames=keys)
        writer.writeheader()
        writer.writerows(all_devs)
        output_file.close()

    # flat_dict['devs'] = pd.DataFrame(devs_list)
    # flat_dict['order_meta'] = pd.DataFrame(map(
    #     lambda x, y, i: {'order': i, 'order_median': x, 'percentile': y[0], 'threshold': y[1]},
    #     emp['order_medians'],
    #     emp['percentiles_and_thresholds'],
    #     range(len(emp['order_medians']))
    # ))
    # d[name] = flat_dict



#  This needs to be massaged into a workable format. The above code is a good scaffold on which to do this.
