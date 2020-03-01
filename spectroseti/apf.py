# apf.py
#
# Nate Tellis 2017
#
#
# Extends class definitions in spectra.py for the High Resolution Echelle
#    spectrometer at the WM Keck observatory.

__author__ = 'nate'
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from pathos.multiprocessing import ProcessingPool as Pool
from tqdm import tqdm

import apfdefinitions as apfdefs
import definitions as defs
import spectra as spec
import utilities as utilities


class APFRedObs(spec.ReducedObs):
    raw_red_correspondence = None

    def __init__(self, run, obs, atlas=spec.WavelengthAtlas()):
        self.devs_set = 0
        self.run = run
        self.obs = obs
        self.loadobservation(run, obs)
        self.atlas = atlas
        self.order_medians = map(np.median, self.counts)
        self.percentiles_and_thresholds = None

    def loadobservation(self, run, obs, deblaze=0):
        """
        Load a reduced APF spectrum into object params
        :param run: 3 character string, i.e. 'aii', denoting observing run
        :param obs: int, denoting observation within a run
        :return: none
        """

        self.devs_set = 0
        self.run = run
        self.obs = obs

        try:
            self.dat = fits.open(apfdefs.spectra_dir_reduced + 'r%(run)s.%(obs)s.fits' % locals())
            self.wav = fits.open(defs.project_root_dir + apfdefs.apf_wavs)
            self.wavs = self.wav[0].data
            self.counts = self.dat[0].data
            if deblaze:
                self.deblaze_orders('savitzky')

        except IOError:
            self.dat = None
            print('Fatal Error - r%(run)s.%(obs)s.fits not found!' % locals())
    def plot_order(self,order,scaling='auto'):

        if scaling == 'auto':
            run=self.run
            obs=self.obs
            plt.figure()
            plt.plot(self.wavs[order,:], self.counts[order,:-1], zorder=1)
            plt.scatter(self.wavs[order,:], self.counts[order,:-1],c='r',s=4., zorder=2)
            if self.percentiles_and_thresholds:
                plt.hlines(self.percentiles_and_thresholds[order][0],self.wavs[order,0],self.wavs[order,-1] , color='k', zorder=3)
                plt.hlines(self.percentiles_and_thresholds[order][1],self.wavs[order,0],self.wavs[order,-1] , color='g', zorder=4)
            plt.title(('Order %(order)s of ' % locals()) + self.dat[0].header['OBJECT'] +
                      ' in r%(run)s.%(obs)s.fits' % locals())
            plt.ylabel('Flux (Counts)')
            plt.xlabel('Wavelength (Angstrom)')

        else:
            raise NotImplementedError

    def plot_all_orders(self):
        for i in range(79):
            run = self.run
            obs = self.obs
            plt.figure(figsize=[15,10])
            plt.plot(self.wavs[i, :], self.counts[i, :-1])
            plt.scatter(self.wavs[i, :], self.counts[i, :-1], c='r', s=4., zorder=2)
            if self.percentiles_and_thresholds:
                plt.hlines(self.percentiles_and_thresholds[i][0],self.wavs[i,0],self.wavs[i,-1] , color='k', zorder=3)
                plt.hlines(self.percentiles_and_thresholds[i][1],self.wavs[i,0],self.wavs[i,-1] , color='g', zorder=4)
            plt.title(('Order %(i)s of ' % locals()) + self.dat[0].header['OBJECT'] +
                      ' in r%(run)s.%(obs)s.fits' % locals())
            plt.ylabel('Flux (Counts)')
            plt.xlabel('Wavelength (Angstrom)')
            plt.savefig(apfdefs.reduced_orders_dir + ('Order %(i)s of ' % locals()) + self.dat[0].header['OBJECT'] + '.png')
            plt.cla()
            plt.close()

    def deblaze_orders(self, method='meanshift', percentile_kernel=501, savitzky_kernel=2001,
                       savitzky_degree=4, perc=50, bstar_correction=None, multi_cores=1):
        if method == 'savitzky':
            deblaze = lambda x: utilities.deblaze(x, method='savitzky', percentile_kernel=percentile_kernel,
                                                  savitzky_kernel=savitzky_kernel,
                                                  savitzky_degree=savitzky_degree, perc=perc)
            self.counts = np.apply_along_axis(deblaze, 1, self.counts)
        elif method == 'percentile':
            deblaze = lambda x: utilities.deblaze(x, method='percentile', percentile_kernel=percentile_kernel,
                                                  savitzky_kernel=savitzky_kernel,
                                                  savitzky_degree=savitzky_degree, perc=perc)
            self.counts = np.apply_along_axis(deblaze, 1, self.counts)
        elif method == 'bstar':
            try:
                bstar_correction.shape
            except AttributeError:
                bstar_correction = np.load(defs.project_root_dir + apfdefs.bstar_correction_dir)
            self.counts = self.counts / bstar_correction

        # Meanshift deblaze only really works once another deblaze has been applied
        elif method == 'meanshift':
            deblaze = lambda x: utilities.deblaze(x, method='meanshift', percentile_kernel=percentile_kernel,
                                                  savitzky_kernel=savitzky_kernel,
                                                  savitzky_degree=savitzky_degree, perc=perc)

            # --------------------------- WARNING - MULTIPROCESSING ---------------
            #  Perform a multicore deblaze
            #  Pathos allows for generic functions to be used (uses dill vs pickle)
            if multi_cores > 1:
                p = Pool(multi_cores)
                pool_output = Pool.map(p, deblaze, self.counts)
                self.counts = np.array(pool_output)

                deblaze = lambda x: utilities.deblaze(x, method='percentile', percentile_kernel=percentile_kernel,
                                                      savitzky_kernel=savitzky_kernel,
                                                      savitzky_degree=savitzky_degree, perc=perc)
                self.counts = np.apply_along_axis(deblaze, 1, self.counts)
            else:
                self.counts = np.apply_along_axis(deblaze, 1, self.counts)

                deblaze = lambda x: utilities.deblaze(x, method='percentile', percentile_kernel=percentile_kernel,
                                                      savitzky_kernel=savitzky_kernel,
                                                      savitzky_degree=savitzky_degree, perc=perc)
                self.counts = np.apply_along_axis(deblaze, 1, self.counts)
        else:
            raise KeyError(
                'The deblaze method you have passed is not implemented. Please pick from savitzky, bstar, and meanshift')

    def loaddevs(self, method='simple', n_mads=5, percentile=75, multi_cores=1):
        if method == 'simple':
            self.devs, self.percentiles_and_thresholds = findhigher(self, n_mads, percentile, atlas=self.atlas)
            self.devs_set = 1
        # An even simpler method that is only 4 pix >1.3*med
        if method == 'multiprocess':
            self.devs, self.percentiles_and_thresholds = findhigher(self, n_mads, percentile, atlas=self.atlas,
                                                                    method='multiprocess',multi_cores=multi_cores)
            self.devs_set = 1
        if method == 'simpler':
            self.devs, self.percentiles_and_thresholds = findhigher(self, n_mads, percentile, atlas=self.atlas,
                                                                    method='basic')
            self.devs_set = 1
        if method == 'MAD':
            raise NotImplementedError

    def red2raw(self, ord, pix):

        try:
            if not self.raw_red_correspondence:
                self.raw_red_correspondence = apfdefs.correspondence
        except IOError:
            print('Fatal Error - apf order correspondence mask not found!')

        return self.raw_red_correspondence[ord, pix]

    #   Warning - this needs a linked raw image
    def dev_has_cr(self, devnum, raw, method='ratio'):
        if self.devs_set:
            ord = self.devs[devnum][0]
            mid = self.devs[devnum][2]+self.devs[devnum][1]/2
            reject = raw.cr_reject(ord, mid,method=method)
            return reject
        else:
            print("Devs have not been set!")

    #   New (DEC 2019) CR rejector
    def reduced_image_cr_check(self, ord, central_pixel, continuum):
        maxpix = np.argmax(self.counts[ord][central_pixel-5:central_pixel+5])
        subset = self.counts[ord][central_pixel+(maxpix-5)-10:central_pixel+(maxpix-5)+10] - continuum
        peak = subset[10]
        subset_norm = subset/peak
        pair_one = subset_norm[9]+subset_norm[11]
        pair_two = subset_norm[8]+subset_norm[12]
        pair_three = subset_norm[7]+subset_norm[13]
        pair_four = subset_norm[6]+subset_norm[14]
        pair_five = subset_norm[5]+subset_norm[15]
        if pair_one > 2.0:
            return True, "pair_one > 2.0"
        elif pair_one < 1.0:
            return True, "pair_one < 1.0"
        elif pair_two < 0.4:
            return True, "pair_two < 0.4"
        elif pair_one < pair_two:
            return True, "pair_one < pair_two"
        elif pair_three < 0.05:
            return True, "pair_three < 0.05"
        elif pair_two < pair_three:
            return True, "pair_two < pair_three"
        elif pair_three < pair_four:
            return True, "pair_three < pair_four"
        elif pair_three < pair_five:
            return True, "pair_three < pair_five"
        else:
            return False, ""



class APFRawObs(spec.RawObs):
    """
    Extends RawObs class for the specifics of the APF
    """
    data = None
    raw_red_correspondence = None
    rawdims = (4608, 2080)

    def __init__(self, run, obs):
        self.load(run, obs)
        self.raw_red_correspondence = np.load(defs.project_root_dir + apfdefs.correspondence)

    def load(self, run, obs):
        """Load a raw observation."""
        raw_fits = fits.open(apfdefs.spectra_dir_raw + 'ucb-%(run)s%(obs)s.fits' % locals())
        self.data = np.rot90(raw_fits[0].data, 3)  # rotation to get raw in line with prefered viewing orientation
        raw_fits.close()

    def show_raw(self):
        fig, ax = plt.subplots()
        ax.imshow(np.transpose(raw[0].data), vmin=np.min(raw[0].data), vmax=np.percentile(raw[0].data, 85))

    def red2raw(self, ord, pix):
        '''try:
            if self.raw_red_correspondence == None:
                self.raw_red_correspondence = np.load(defs.project_root_dir + apfdefs.correspondence)
        except IOError:
            print('Fatal Error - apf order correspondence mask not found!')'''

        return self.raw_red_correspondence[ord, pix]

    def retrieve_subset(self, ord, index, yradius=5, xradius=5):
        """
        Retrieves a rectangular subset of the raw image at
        the stellar ridge corresponding to some order/pixel pair
        """
        central = self.red2raw(ord, index)
        ylow = 0 if central - yradius < 0 else central - yradius
        yhigh = self.rawdims[1] if central + yradius > self.rawdims[1] else central + yradius
        xlow = 0 if index - xradius < 0 else index - xradius
        xhigh = self.rawdims[0] if index - xradius > self.rawdims[0] else index + xradius
        return self.data[int(ylow):int(yhigh), int(xlow):int(xhigh)]

    def cr_reject(self, order, pix, method='ratio',xradius=20,yradius=12):
        if method == 'std':
            postage_stamp = self.retrieve_subset(order, pix, yradius=yradius, xradius=xradius)
            means = np.repeat(np.reshape(np.apply_along_axis(
                lambda x: np.mean(x[x<np.percentile(x,90)]),1,postage_stamp),
                (yradius*2,1)),repeats=xradius*2,axis=1)
            #
            # plt.imshow(means, cmap='viridis')
            # plt.colorbar()
            stds = np.repeat(np.reshape(np.apply_along_axis(
                lambda x: np.std(x[x<np.percentile(x,90)]),1,postage_stamp),
                (yradius*2,1)),repeats=xradius*2,axis=1)
            #
            # plt.imshow(stds, cmap='viridis')
            # plt.colorbar()
            comp = (postage_stamp - means) / stds
            #plt.imshow(comp, cmap='viridis')
            #plt.colorbar()
            #print(comp[:,xradius-3:xradius+3])
            return np.max(comp[:,xradius-3:xradius+3])
            # if np.any(comp[:,xradius-3:xradius+3] > 12):
            #     return
            # else:
            #     return 'Does not appear to be a cosmic ray.'
        elif method == 'ratio':
            postage_stamp = self.retrieve_subset(order, pix, yradius=yradius, xradius=12)

            one_step,two_step, maxval = utilities.compute_maxpix_deviance(postage_stamp)
            max1spatial = np.sum(np.array(one_step[:2]) < 0.20)
            max2spatial = np.sum(np.array(two_step[:2]) < 0.06)
            max1spec = np.sum(np.array(one_step[2:]) < 0.20)
            max2spec = np.sum(np.array(two_step[2:]) < 0.06)
            compval = max1spatial + max2spatial + max1spec +max2spec
            return compval

    def run_obs_to_filename(self, run, obs):
        """Simple utility to translate from run/obs int pair to raw filename."""
        return 'ucb-%(run)s%(obs)s.fits' % locals()


# DEPRECATED
'''def populateapfobsdatabase(filename='apf_obs.db'):
    filenames = listdir(apfdefs.logsheet_dir)
    log_accumulator = []
    for log in filenames: # Replace with all logsheets
        run = log.split('.')[0]
        f = open(apfdefs.logsheet_dir+log)
        logs = f.readlines()
        for i in range(7,13):
            try:
                if logs[i][:3].isdigit():
                    break
            except IndexError:
                continue
        log_accumulator += list(map(lambda x: str.split(run+' '+x,maxsplit=12), logs[i:]))
        f.close()
        length = len(sorted(log_accumulator,key=len, reverse=True)[0])
        acc_array = np.array([xi+[None]*(length-len(xi)) for xi in log_accumulator])'''


# Turned off AltMinThresh for now. See what this does.
def find_deviations(ords, wavs, order, perc=75, n_mads=5, alt_min_thresh=0, atlas=spec.WavelengthAtlas(), out=[],
                    npix=3, acc=[], precomputed_percentile = None, precomputed_threshold = None):
    """
        The meat of the laser search pipeline, this function finds all pixel regions that deviate
        over an entire spectroscopic order.

    :param ords: Order counts , numpy.ndarray
    :param wavs: Corresponding wavelength range, numpy.ndarray
    :param order: Which order to find deviations in?
    :param simulator: DEPRECATED - for adding simulated lasers to a run
    :param perc: What(th) percentile should the MAD be computed above? For continuum determination
    :param n_mads: How many MADs should the code look for deviations above?
    :param atlas: Wavelength atlas for night sky and emission line removal
    :param out: Can pass a list for accumulation
    :param npix: How many pixels to demand consecutive
    :return: list of lists with format [[order,start,len,nth_percentile, threshold set, mean deviate pixel
                                        median deviate pixel, midpt pixel value], ...]
    """
    l = len(ords[0])
    if not precomputed_percentile:
        percentile = utilities.getpercentile(ords[order], perc)
    else:
        percentile = precomputed_percentile
    if not precomputed_threshold:
        threshold = utilities.findthresh(ords[order] - percentile)
    else:
        threshold = precomputed_threshold
    th2 = 100 + percentile  # TODO - fix second thresholding
    if percentile < 100 and threshold < th2:
        final_threshold = percentile + n_mads * threshold * 1.5
        secondary_threshold = percentile + n_mads * threshold / 2.0
    else:
        final_threshold = percentile + n_mads * threshold
        secondary_threshold = percentile + n_mads * threshold / 2.0
    # Experimental basic thresholding so that MADS is not far too high
    if percentile > 500 and n_mads * threshold > 0.3 * percentile and alt_min_thresh:
        final_threshold = 1.3 * percentile
        secondary_threshold = 1.15 * percentile
    # Accumulator for testing whole-dataset thresholding
    acc.append([percentile, final_threshold])
    contig = utilities.finddeviates(ords[order], final_threshold, npix=npix)
    if len(contig) != 0:
        for x in contig:
            midpt = x[0] // 2 + x[1]
            if x[1] > 500 and x[0] + x[1] < 4108:# \ COMMENTED OUT THE IGNORED WAVELENGTHS
                # and not hires_ignored_wavelengths(wavs[order][midpt]) \
                # and not list(atlas.ns_lookup(wavs[order][midpt])) and not has_singularity(ords[order]):
                deviation_pixels = ords[order][x[1]:x[1] + x[0]]
                out.append(
                    [order, x[0], x[1], float(percentile), float(threshold), float(np.mean(deviation_pixels)),
                     float(np.median(deviation_pixels)), float(wavs[order][midpt]), max(deviation_pixels)])
    return out


def find_deviations_basic(ords, wavs, order, perc=75, n_mads=5, alt_min_thresh=1, atlas=spec.WavelengthAtlas(), out=[],
                          npix=3, acc=[]):
    """
        The meat of the laser search pipeline, this function finds all pixel regions that deviate
        over an entire spectroscopic order.

    :param ords: Order counts , numpy.ndarray
    :param wavs: Corresponding wavelength range, numpy.ndarray
    :param order: Which order to find deviations in?
    :param simulator: DEPRECATED - for adding simulated lasers to a run
    :param perc: What(th) percentile should the MAD be computed above? For continuum determination
    :param n_mads: How many MADs should the code look for deviations above?
    :param atlas: Wavelength atlas for night sky and emission line removal
    :param out: Can pass a list for accumulation
    :param npix: How many pixels to demand consecutive
    :return: list of lists with format [[order,start,len,nth_percentile, threshold set, mean deviate pixel
                                        median deviate pixel, midpt pixel value], ...]
    """
    l = len(ords[0])
    percentile = utilities.getpercentile(ords[order], perc)
    threshold = utilities.findthresh(ords[order] - percentile)
    final_threshold = 1.5 * percentile
    # Accumulator for testing whole-dataset thresholding
    acc.append([percentile, final_threshold])
    contig = utilities.finddeviates(ords[order], final_threshold, npix=npix)
    if len(contig) != 0:
        for x in contig:
            midpt = x[0] // 2 + x[1]
            if x[1] > 500 and x[0] + x[1] < 4108 and x[0] < 10:  # \ COMMENTED OUT THE IGNORED WAVELENGTHS
                # and not hires_ignored_wavelengths(wavs[order][midpt]) \
                # and not list(atlas.ns_lookup(wavs[order][midpt])) and not has_singularity(ords[order]):
                deviation_pixels = ords[order][x[1]:x[1] + x[0]]
                out.append(
                    [order, x[0], x[1], float(percentile), float(threshold), float(np.mean(deviation_pixels)),
                     float(np.median(deviation_pixels)), float(wavs[order][midpt]), max(deviation_pixels)])
    return out


def findhigher(obs, n_mads, perc, atlas=spec.WavelengthAtlas(), method='original', multi_cores = 1):
    """
    Finds potential laser candidates using Median
    Absolute Deviation method in find_deviations,
    potentially adding simulated signals as well.

    First bins spectrum, having subtracted 75th percentile
    then finds median 4-pix deviation, uses this as threshold
    by which to search for every multi-pixel deviation
    of at least 4 pixels


    :param atlas: Wavelength atlas for throwing out night sky line frequencies
    :param obs: HiresObs instance
    :param n_mads: Number of MADs above %ile to use as threshold
    :param perc: Percentile to use for stellar continuum
    :return: list of lists with format [[order,start,len,nth_percentile, threshold set, mean deviate pixel
                                        median deviate pixel, midpt pixel value], ...]
    """
    out = []
    per_thr_accumulator = []
    cts = obs.counts
    wavs = obs.wavs
    if method == 'original':
        for i in tqdm(range(79), leave=False, miniters=8):
            out = find_deviations(cts, wavs, i, perc=perc, n_mads=n_mads,
                                  atlas=atlas, out=out, acc=per_thr_accumulator)
        return out, per_thr_accumulator


    # --------------------------- WARNING - MULTIPROCESSING ---------------
    #  Perform a multicore search for laser lines by pulling out the percentile, threshold computation
    #  Pathos allows for generic functions to be used (uses dill vs pickle)
    #  Take care not to overwhelm shared resources

    elif method == 'multiprocess':
        # Lambda function denoting the xth percentile (continuum value)
        compute_percentile = lambda x: utilities.getpercentile(cts[x],perc)
        p = Pool(multi_cores)
        pool_output = Pool.map(p, compute_percentile, range(79))
        percentiles = np.array(pool_output)
        for i in tqdm(range(79), leave=False, miniters=8):
            out = find_deviations(cts, wavs, i, perc=perc, n_mads=n_mads,
                                        atlas=atlas, out=out, acc=per_thr_accumulator,
                                        precomputed_percentile=percentiles[i])
        return out, per_thr_accumulator
    # Deprecated, used for a simplistic > 1.3* continuum threshold
    elif method == 'basic':
        for i in tqdm(range(79), leave=False, miniters=8):
            out = find_deviations_basic(cts, wavs, i, perc=perc, n_mads=n_mads,
                                        atlas=atlas, out=out, acc=per_thr_accumulator, npix=4)
        return out, per_thr_accumulator
    else:
        raise NameError('Incorrect method keyword')
