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
import os
import shutil
import pandas as pd
from xarray import Dataset
import glob


candidates = pd.read_csv('../data/allAPFCandidateMetadataOct04_filtered.csv')
all_image_loc = '/media/nate/DATA/code/output/las_images/allCandidateImages/'
dest_dir = '/media/nate/DATA/code/oustput/las_images/intenseCandidateImages/'
minwav = 3735.

candidates = candidates.loc[~candidates['target_name'].isin(['NarrowFlat', 'test', 'Dark'])]
# loop over the bins ( min wav + 5 A increments)
for i in range(1295):
    subset = candidates.loc[(candidates['central_wav'] > (minwav+i*5.)) & (candidates['central_wav'] <= (minwav+(i+1)*5.))]
    subset_nocr = subset.loc[subset['cosmic_reject_value']<2]
    filenames = subset_nocr.sort_values('intensity', ascending=False)['img_prefix'][0:3].values
    print filenames
    for counter, filename in enumerate(filenames):
        for file in glob.glob(all_image_loc+r'%(filename)s_*.jpg' % locals()):
            print(file)
            printable_wav = minwav+(i+0.5)*5.
            shutil.copy(file, dest_dir+ '%(printable_wav)s_%(counter)s_' %locals()+filename+'.jpg')


    #sort each bin by intensity

#copy the images to their own folder, wioth wavelength filename preset