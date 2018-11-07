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






# Modulate these
medfilt_kernel_size = 501
savitzky_kernel_size = 51
savitzky_order = 4

mdf = lambda x: sg.medfilt(x,kernel_size=medfilt_kernel_size)
bs_medfilt = np.apply_along_axis(mdf, 1, bstar_mean)

bs_medfilt_savitzky = np.apply_along_axis(lambda x: util.savitzky_golay(x, savitzky_kernel_size, savitzky_order),
                                          1, bs_medfilt)



#Use this corrector for ords up to 10?
inds = np.array([8,9,12,13,14,15])



#iterate over:


for i in range(15,70):
    inds = np.arange(i-4,i+5)





for i in range(77):
    order = i
    testorder = red.counts[order,:-1]
    testwavs = red.wavs[order,:]





    # Works very well for order 20
    # corrector_ord7 =(np.max(testorder[2200:-2200]*2)/(1+np.power((testwavs-testwavs[int(len(testwavs)/2.15)])/38,2))) - np.max(testorder[2200:-2200])

    # Works very well for order 7
    # corrector_ord7 =((testwavs[2400])/(1+np.power((testwavs-testwavs[len(testwavs)/2])/22,2))-1100) / 1.25


    # Order 6
    corrector_ord7 =((testwavs[2400])/(1+np.power((testwavs-testwavs[int(len(testwavs)/2.1)])/22,2))-1100) / 1.35

    secondary_corrector = bs_medfilt_savitzky[order,:-1]

    # For order 7 (CA2 H
    # secondary_corrector[100:2500] = corrector_ord7[100:2500] / np.percentile(testorder,98)*np.max(bs_medfilt[order,:-1])



    # For order 6 (H and K)
    secondary_corrector = corrector_ord7 / np.percentile(testorder,98)*np.max(bs_medfilt[order,:-1])


    plt.figure(1)
    plt.plot(testwavs,testorder)
    plt.plot(testwavs, bs_medfilt_savitzky[order,:-1]/np.max(bs_medfilt_savitzky[order,:-1])* np.percentile(testorder,98))
    plt.plot(testwavs, bs_medfilt[order,:-1]/np.max(bs_medfilt[order,:-1])* np.percentile(testorder,98))
    # plt.plot(testwavs,corrector_ord7)
    plt.scatter(testwavs,testorder,c='r',s=4.)


    # plt.figure(2)
    # plt.plot(testwavs,testorder/secondary_corrector)
    # # plt.scatter(testwavs,testorder/bs_medfilt_savitzky[order,:-1],c='r',s=4.)
    # plt.scatter(testwavs,testorder/secondary_corrector,c='r',s=4.)

    plt.savefig("bstar_order_%(i)s.png" % locals())
    plt.cla()