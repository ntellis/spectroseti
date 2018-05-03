# injector.py
#
# Nate Tellis 2017
#
#
# Contains classes and methods for injectiona nd recovery testing
#
#

import spectroseti.output

__author__ = 'nate'


import definitions as defs
import apf as apf
import apfdefinitions as apfdefs
#import output
import numpy as np
from tqdm import tqdm
import random


gaussian_kernel_4 = np.array([0.02169646, 0.08608803, 0.2514724, 0.54092311, 0.85695116,
                                  1., 0.85957485, 0.54424083, 0.25379, 0.0871478, 0.02203095])

class Injector:

    fwhm = 4.2 # CHECK THIS
    injected_lasers = None
    injected = False
    spectrum = None
    spectrum_backup = None


    def __init__(self):
        pass


    def setparams(self,fwhm=4.2,):
        self.fwhm=fwhm

    def set_spectrum(self,spectrum):
        self.spectrum = spectrum
        #deep copy spectrum for backup
        #Typecheck spectrum here


    # Generates flux levels and locations, stores in injected_lasers
    def gen_lasers(self,per_order = 10, fluxes = [], autodetermine_flux = True):
        if autodetermine_flux:
            if self.spectrum:
                # commpute flux levels from spectrum
                # base this on order medians?
                pass
            else:
                fluxes = [10, 20, 50, 100, 500, 1000, 5000]
        acc=[]
        #this could be cleaner
        for i in range(79):
            for j in range(per_order):
                # Maybe put a flux in here too?

                acc.append([i,random.randint(500,4300), random.choice(fluxes)])
        self.injected_lasers = acc





    def inject_lasers(self):
        #check if the lasers have been generated
        if self.injected_lasers and self.spectrum:
            for laser in self.injected_lasers:
                kernel = gen_laser_profile(power=laser[2])
                #Update this so that the backup remains
                low = laser[1]-6
                high = laser[1]+6
                print(self.spectrum.counts[laser[0],low:high])
                self.spectrum.counts[laser[0],low:high] = \
                    self.spectrum.counts[laser[0],low:high] + kernel
                print(self.spectrum.counts[laser[0], low:high])
                print('Next')
            self.injected=True
        else:
            print('Either spectrum is not set or lasers are not chosen')

    # Test for recovery of simulated signals
    def recovery(self, deblaze='both', bypass_threshold=True, repetitions=100, per_order=10,):
        pass

        # This should take deblazing parameters, repetitions

        # Should run the target without any injection, collect positives, mask out these locations
        # SHould also mask out the locations of stellar emission lines
        # Should also mask out night sky lines?

        # Run, record
            # Inject, record injections
                # Run, compare recovery, inc overlapping injections

        # log to a file or to terminal

        # How to speed this up? It will be very hard to run thousands of trials at a minute each
        # One possibility is to inject the lasers at the end, pretending that it didn't affect deblazing?
        # This would allow us to run a bunch of lasers after deblazing, and threshold have been determined.
        # in fact it makes more sense to add lasers at this stage (post - deblazing) due to flattening
        '''
        So the flow is:
        deblaze savitzky
        opt: store (deep copy) spectrum
        deblaze meanshift
        run search
        opt: run search on savitzky spectrum
        '''

def gen_laser_profile(power=1.0):
    kernel = np.zeros(12)
    pixel_offset_weight = random.random()
    # Assign most of the counts to one pixel, the rest next door
    kernel[:-1] = kernel[:-1] + gaussian_kernel_4 * power * pixel_offset_weight
    kernel[1:] = kernel[1:] + gaussian_kernel_4 * power * (1.0 - pixel_offset_weight)
    return kernel