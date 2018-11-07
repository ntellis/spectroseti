# Bstar_corrections.py
# Nathaniel Tellis, Jan 2018
#
# Generates and saves a .npy that contains apf corrections of three (five?) averaged b-star continua for deblazing


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import seaborn as sb
import spectroseti.apf as apf
import spectroseti.utilities as util

import scipy.signal as sg

from scipy.interpolate import CubicSpline
from sklearn.neighbors import KernelDensity
from sklearn.cluster import MeanShift, estimate_bandwidth

bstar1 = fits.open('/media/nate/DATA/Spectra/apf_bstar/ralh.272.fits')[0].data
bstar2 = fits.open('/media/nate/DATA/Spectra/apf_bstar/rayf.213.fits')[0].data
bstar3 = fits.open('/media/nate/DATA/Spectra/apf_bstar/rayi.237.fits')[0].data

# create a 4096x79x3 array, and take (median? mean?

bstar1 = np.apply_along_axis(util.median_of_one, 1, bstar1)
bstar2 = np.apply_along_axis(util.median_of_one, 1, bstar2)
bstar3 = np.apply_along_axis(util.median_of_one, 1, bstar3)

bstar_mean = (bstar1 + bstar2 + bstar3) / 3.

bs_medfilt = np.apply_along_axis(sg.medfilt(kernel_size=501), 1, bstar_mean)

bs_medfilt_savitzky = np.apply_along_axis(lambda x: util.savitzky_golay(x, 51, 4), 1, bs_medfilt)
