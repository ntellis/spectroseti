import spectroseti.apf as apf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import seaborn as sb
import spectroseti.apf as apf
import spectroseti.utilities as util

import scipy.signal as sg



red = apf.APFRedObs('bac', 249)


bstar1 = fits.open('/media/nate/DATA/Spectra/apf_bstar/ralh.272.fits')[0].data
bstar2 = fits.open('/media/nate/DATA/Spectra/apf_bstar/rayf.213.fits')[0].data
bstar3 = fits.open('/media/nate/DATA/Spectra/apf_bstar/rayi.237.fits')[0].data

# create a 4096x79x3 array, and take (median? mean?

bstar1 = np.apply_along_axis(util.median_of_one, 1, bstar1)
bstar2 = np.apply_along_axis(util.median_of_one, 1, bstar2)
bstar3 = np.apply_along_axis(util.median_of_one, 1, bstar3)

bstar_mean = (bstar1 + bstar2 + bstar3) / 3.

bstar_mean_medianfiltered = np.array(bstar_mean)
bstar_mean_meanfiltered = np.array(bstar_mean)

for i in range(77):
    omit_ords = np.array([0,1,2,3,4,5,6,7,10,11,16,17,28,34,45,46,50,52,53,56,63,64,65,66])
    # Here we pick the corrector to use
    if i < 13:
        inds = np.array([8, 9, 12, 13, 14, 15])
    elif i >= 59:
        inds = np.array([54,55,57,58,59,60,61,62])
    else:
        inds = np.arange(i - 4, i + 5)

    ords_to_use = np.setdiff1d(inds,omit_ords)

    bstar_mean_medianfiltered[i,:] = np.median(bstar_mean[ords_to_use],axis=0)
    bstar_mean_meanfiltered[i,:] = np.mean(bstar_mean[ords_to_use],axis=0)


# Modulate these

medfilt_kernel_size = 501
savitzky_kernel_size = 51
savitzky_order = 4

mdf = lambda x: sg.medfilt(x, kernel_size=medfilt_kernel_size)
bs_medfilt = np.apply_along_axis(mdf, 1, bstar_mean_medianfiltered)
bs_m_medfilt = np.apply_along_axis(mdf, 1, bstar_mean_meanfiltered)
bs_medfilt_savitzky = np.apply_along_axis(lambda x: util.savitzky_golay(x, savitzky_kernel_size, savitzky_order),
                                          1, bs_medfilt)


bs_medfilt_savitzky = np.array(bs_medfilt_savitzky)
bs_medfilt_savitzky[:,:-20] = np.array(bs_medfilt_savitzky[:,20:])


bstar_correction = bs_medfilt_savitzky
np.save("bstar_correction.npy",bstar_correction )

for i in range(77):
    order = i

    testorder = red.counts[order, :-1]
    testwavs = red.wavs[order, :]




    plt.figure(1)
    plt.plot(testwavs,testorder)
    plt.plot(testwavs, bs_medfilt_savitzky[order,:-1]/np.max(bs_medfilt_savitzky[order,:-1])* np.percentile(testorder,98))
    plt.plot(testwavs, bs_medfilt[order,:-1]/np.max(bs_medfilt[order,:-1])* np.percentile(testorder,98))
    # plt.plot(testwavs,corrector_ord7)
    plt.scatter(testwavs,testorder,c='r',s=4.)
    plt.savefig("bstar_order_%(i)s.png" % locals())
    plt.cla()
    if i<=100:
        plt.figure(2,figsize=[10,7])
        plt.plot(testwavs,testorder[:]/bs_medfilt_savitzky[order,:-1])
        plt.scatter(testwavs,testorder[:]/bs_medfilt_savitzky[order,:-1],c='r',s=4.)
        plt.savefig("bstar_order_%(i)s_corrected.png" % locals())
        plt.cla()
    # elif i>30 and i<53:
    #
    #     plt.figure(2, figsize=[10, 7])
    #     plt.plot(testwavs[:-10], testorder[:-10] / bs_medfilt_savitzky[order, 10:-1])
    #     plt.scatter(testwavs[:-10], testorder[:-10] / bs_medfilt_savitzky[order, 10:-1], c='r', s=4.)
    #     plt.savefig("bstar_order_%(i)s_corrected.png" % locals())
    #     plt.cla()
    # else:
    #     plt.figure(2, figsize=[10, 7])
    #     plt.plot(testwavs[:], testorder[:] / bs_medfilt_savitzky[order, :-1])
    #     plt.scatter(testwavs[:], testorder[:] / bs_medfilt_savitzky[order, :-1], c='r', s=4.)
    #     plt.savefig("bstar_order_%(i)s_corrected.png" % locals())
    #     plt.cla()