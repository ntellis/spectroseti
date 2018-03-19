# Basic script for accumulating run statistics


import matplotlib
import os
matplotlib.use('TkAgg')
import multiprocessing as mp

import spectroseti.runner as runner


spec = os.listdir('/media/nate/DATA/Spectra/apf')
# just grabs the 'bac' run
# TODO - turn this into a "search by run" method in LaserSearch
observations = map(lambda x: [x[0][1:], int(x[1])],
                    map(lambda x: x.split('.')[:2],
                        filter(lambda x: x.__contains__('bac'), spec)))

observations.sort()

LS = runner.LaserSearch()

LS.search_multiple(observations[24:],output_pngs=0,multi_cores=mp.cpu_count(), stats=1)