# apf.py
#
# Nate Tellis 2017
#
#
# Extends class definitions in spectra.py for the High Resolution Echelle
#    spectrometer at the WM Keck observatory.
import numpy as np
from matplotlib import pyplot as plt

__author__ = 'nate'


def view_dev(spec,devnum=0,raw=None, save='', nmads=5):
    # if Raw is supplied, plots raw postage stamp

    # This prevents the figure from popping up on screen before save

    import apfdefinitions
    ord = spec.devs[devnum][0]
    try:
        plt.close()
        r = spec.run
        o = spec.obs
        central = spec.devs[devnum][2] + spec.devs[devnum][1]//2
        low = spec.devs[devnum][2] + spec.devs[devnum][1]//2 - 48
        high = spec.devs[devnum][2] + spec.devs[devnum][1]//2 + 73
        halfwidth = spec.devs[devnum][1]//2 + 1
        mid = spec.devs[devnum][2] + spec.devs[devnum][1]//2
        centralwav = spec.wavs[ord, mid]
        plt.figure(figsize=(16, 11), dpi=80)
        plt.title('Deviation number %(devnum)s in r%(r)s.%(o)s.fits, target: . Order %(ord)s, central pixel %(central)s' % locals())
        ax1 = plt.subplot(311)
        ax1.plot(spec.wavs[ord],spec.counts[ord,:-1])
        ax1.plot(spec.wavs[ord, mid - halfwidth:mid + halfwidth+1], spec.counts[ord, mid - halfwidth:mid + halfwidth+1], color='r')
        plt.hlines(spec.devs[devnum][3], spec.wavs[ord, 0], spec.wavs[ord, -2], colors='r',
                   label='meanshift result')
        plt.hlines(spec.devs[devnum][3] + spec.devs[devnum][4] * nmads, spec.wavs[0, low], spec.wavs[ord, -2],
                   colors='g', label='Threshold (ms+5*MAD)')
        ymax = np.max([np.percentile(spec.counts[ord,:-1],99),np.max(spec.counts[ord,mid-10:mid+10])*1.05])
        ax1.set_ylim([0,ymax])
        ax1.set_xlim(spec.wavs[ord][0]-.5,spec.wavs[ord][-1]+.5)
        plt.xlabel('Wavelength (Ang)')
        plt.ylabel('Counts per pixel')
        ax2 = plt.subplot(312)
        ax2.plot(spec.wavs[ord, low:high], spec.counts[ord, low:high])
        ax2.scatter(spec.wavs[ord, low:high], spec.counts[ord, low:high], color='b')
        ax2.plot(spec.wavs[ord, mid - halfwidth:mid + halfwidth + 1],
                 spec.counts[ord, mid - halfwidth:mid + halfwidth + 1], color='r')
        ymax = np.max([np.percentile(spec.counts[ord, :-1], 99), np.max(spec.counts[ord, mid - 10:mid + 10]) * 1.25])
        ax2.set_ylim([0, ymax])
        plt.hlines(spec.devs[devnum][3], spec.wavs[ord, low], spec.wavs[ord, high], colors='r',
                   label='meanshift result')
        plt.hlines(spec.devs[devnum][3] + spec.devs[devnum][4] * nmads, spec.wavs[ord, low], spec.wavs[ord, high],
                   colors='g', label='Threshold (ms+5*MAD)')
        ax2.set_xlim(spec.wavs[ord][low] - .05, spec.wavs[ord][high] + .05)
        plt.xlabel('Wavelength (Ang)')
        plt.ylabel('Counts per pixel')
        reject = 'no_raw'
        if raw:
            reject = raw.cr_reject(ord, mid)
        ax = plt.gca()
        targname = spec.dat[0].header['OBJECT']
        ax.set_title(
            'Deviation number %(devnum)s in r%(r)s.%(o)s.fits, target: %(targname)s   Order %(ord)s, central pixel %(central)s, reject value: %(reject)s' % locals())
        if raw:
            plt.subplot(313)
            subplot_xrad = 50
            subplot_yrad = 15
            subset = raw.retrieve_subset(ord,mid,yradius=subplot_yrad,xradius=subplot_xrad)
            subsub = raw.retrieve_subset(ord,mid,yradius=4,xradius=6)
            ax3=plt.imshow(subset,interpolation='nearest', vmax = np.max(subsub),
                           vmin = np.min(subset),cmap='plasma',
                           extent=[spec.wavs[ord,mid-subplot_xrad],spec.wavs[ord,mid+subplot_xrad],
                                   0-subplot_yrad//2,subplot_yrad//2], aspect='auto')
            plt.colorbar()
        plt.tight_layout()
        if save:
            plt.savefig(save+'r%(r)s.%(o)s_dev%(devnum)s_%(centralwav)s.png' % locals())
            plt.close()
    except IndexError:
        print('Out of bounds')
    except ValueError:
        print('Value error caught - probably a mismatch size')