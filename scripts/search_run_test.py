import matplotlib
from pathos.multiprocessing import ProcessingPool as Pool

matplotlib.use('TkAgg')

import spectroseti.runner as runner

LS = runner.LaserSearch()

#LS.search_multiple(observations, output_pngs=1, number_mads=8)

results = LS.search_run(run='awx', multi_cores=0, output_pngs=1, db_write=1,number_mads=10, search_title='errorfixing')
