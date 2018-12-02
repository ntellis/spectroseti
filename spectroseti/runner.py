# apf.py
#
# Nate Tellis 2017
#
#
# Extends class definitions in spectra.py for the High Resolution Echelle
#    spectrometer at the WM Keck observatory.
import spectroseti.output

__author__ = 'nate'


import definitions as defs
import apf as apf
import apfdefinitions as apfdefs
#import output
import numpy as np
from tqdm import tqdm
import pandas as pd
from os import listdir, mkdir
from pathos.multiprocessing import ProcessingPool as Pool

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InvalidTargetError(Error):
    """Exception raised when a target that Laser Search does not run on is accessed.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message



# highest level laser search class. Should contain helper methods to extract specific targets etc
class LaserSearch():

    raw_directory = apfdefs.spectra_dir_raw
    reduced_directory = apfdefs.spectra_dir_reduced
    bstar_correction = np.load(defs.project_root_dir + apfdefs.bstar_correction_dir)

    # Default initialization
    # sets raw, red directories to those in apf_def
    def __init__(self, raw_dir = apfdefs.spectra_dir_raw, red_dir = apfdefs.spectra_dir_reduced):
        self.raw_directory = raw_dir
        self.reduced_directory = red_dir
        pass




    def search_one(self, run, observation, load_devs_method='simple',number_mads=5,deblaze_method='percentile',
                   search_percentile=75, multi_cores=1):
        # don't need the raw at first
        # raw = apf.APFRawObs(run, observation)
        reduced_spectrum = apf.APFRedObs(run, observation)

        if (reduced_spectrum.dat[0].header['TOBJECT'] in apfdefs.ignore_targets_apf)\
                or (reduced_spectrum.dat[0].header['OBJECT'] in apfdefs.ignore_targets_apf):
            raise InvalidTargetError(None, reduced_spectrum.dat[0].header['TOBJECT'])

        # Now first deblaze (Savizky works better on lower orders than b-star)
        reduced_spectrum.deblaze_orders(method=deblaze_method,bstar_correction=self.bstar_correction)

        # Make a copy of the (bstar-only) deblazed spectrum
        #bstar_deblazed = copy.deepcopy(reduced_spectrum)
        # print('Beginning Meanshift deblazing')
        # # Meanshift deblaze the reduced spectrum
        # reduced_spectrum.deblaze_orders(method='meanshift', multi_cores=multi_cores)
        # # Load deviations with the meanshift method
        # loaddevs -> findhigher -> find_deviations -> getpercentile (has meanshift method)
        reduced_spectrum.loaddevs(method=load_devs_method,n_mads=number_mads,
                                  percentile=search_percentile, multi_cores=multi_cores)
        # Here we go back and check the bstar spectrum for the same positives
        # One way to proceed is:
        #   compute perc from bstar_deblazed
        #   compute thresh
        #   ensure three pixels are > perc_b+n_mads*thresh_b

        return reduced_spectrum


    def search_multiple(self, observations, output_pngs=0, logfile=0, deblaze_method='percentile',
                        db_write=0, stats=0, multi_cores=1,number_mads = 5,quiet=1, search_title="TitleUnset"):
        # observations expects a tuple of run,obs pairs
        # setup directories, filenames, local accumulator variables etc
        ctr = 1

        if stats:
            deviation_dict = dict()
            perc_thresh_dict = dict()
        if output_pngs or logfile:
            try:
                mkdir(apfdefs.laser_search_run_dir + search_title)
            except OSError:
                pass

        # a little pseudocodey
        for observation in observations:
            print(observation)
            if observations[0][0] =='r':
                fn_split = observation.split('.')
                run = fn_split[0][1:]
                obs = fn_split[1]
            else:
                run = observation[0]
                obs = observation[1]
            try:
                if multi_cores>1:
                    method = 'multiprocess'
                else:
                    method = 'simple'
                reduced_spectrum = self.search_one(run, obs, load_devs_method=method, deblaze_method=deblaze_method,
                                                   multi_cores=multi_cores,number_mads=number_mads)
            # TODO - this does not work for now
            except InvalidTargetError as err:
                print('Attempted to perform search on:    ' + err.message)
                print('Skipping....\n')
                continue
            except Error:
                print("An error occurred!")
                continue

            if output_pngs:
                # generate and save all output
                raw = None
                try:
                    raw = apf.APFRawObs(run,obs)
                except:
                    raw = None
                ndev = len(reduced_spectrum.devs)

                savedir = apfdefs.laser_search_run_dir + search_title + '/'
                print('Writing output images to '+ apfdefs.laser_search_run_dir + search_title)
                np.savetxt(savedir + ('r%(run)s_%(obs)s_percentiles_and_thresholds.txt' % locals()),reduced_spectrum.percentiles_and_thresholds )
                if quiet:
                    for i in range(ndev):
                        ord = reduced_spectrum.devs[i][0]
                        mid = reduced_spectrum.devs[i][2] + reduced_spectrum.devs[i][1] // 2
                        reject = raw.cr_reject(ord, mid)
                        # if reject >=2:
                        #     continue

                        print(reduced_spectrum.devs[i][-1])
                        if (reduced_spectrum.devs[i][-1] > 6868. and reduced_spectrum.devs[i][-1] < 6885):
                            print('REJECTED')
                            continue
                        elif (reduced_spectrum.devs[i][-1] > 7595. and reduced_spectrum.devs[i][-1] < 7618.):
                            print('REJECTED')
                            continue
                        spectroseti.output.view_dev(reduced_spectrum, devnum=i, raw=raw, save=savedir,nmads=number_mads)
                else:
                    for i in tqdm(range(ndev), miniters=int(ndev/10)):
                        spectroseti.output.view_dev(reduced_spectrum, devnum=i, raw=raw, save=savedir)
                pass

            #TODO this is first priority
            if logfile:
                # Accumulate statistics for logfile
                pass

            if db_write:
                # Save to database
                pass

            if stats:
                title = tuple(observation)
                perc_thresh_dict[title] = reduced_spectrum.percentiles_and_thresholds
                deviation_dict[title] = reduced_spectrum.devs
                # Just save every iteration in case it breaks
                np.save(defs.project_root_dir+'/spectroseti/data/Percentiles_and_thresholds2.npy', perc_thresh_dict)
                np.save(defs.project_root_dir+'/spectroseti/data/Candidates2.npy', deviation_dict)




        # write logfile
        if logfile:
            #save logfile to logfile directory, as well as search run directory
            pass

    def search_run(self,run='bac', output_pngs=0, logfile=0,
                        db_write=0, stats=0, multi_cores=1,number_mads = 5, search_title="TitleUnset" ):
        # Look in the reduced dir for files corresponding to this run
        all_reduced = listdir(self.reduced_directory)
        files_to_search = [fn for fn in all_reduced if fn[1:4] == run]
        p = Pool(multi_cores)
        search_multi = lambda x: self.search_multiple([x], output_pngs=output_pngs, logfile=logfile, db_write=db_write,
                                                      stats=stats, number_mads=number_mads, search_title=search_title)

        pool_output = Pool.map(p, search_multi, files_to_search)
        return pool_output




def RunObsFromLogsheet(filename):
    #check fixrows
    #check header
    pd.read_fwf(filename, skiprows=11, header=None,
                names=['Obs', 'Target', 'Iodine', 'UTC', 'Exptime', 'Deck', 'Counts', 'HA', 'AZ/EL', 'SpecTemp',
                       'DomeTemp', 'Seeing'])
    pass

