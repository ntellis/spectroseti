import matplotlib
from pathos.multiprocessing import ProcessingPool as Pool

matplotlib.use('TkAgg')

import spectroseti.runner as runner

observations = [['awx', 222]]  # ,[],['awx',223],['awx',224]]
bacrun = [['bac', 196],
          ['bac', 273]]

LS = runner.LaserSearch()

#LS.search_multiple(observations, output_pngs=1, number_mads=8)

p = Pool(6)
search_multi = lambda x: LS.search_multiple([x], output_pngs=1,number_mads=10)

pool_output = Pool.map(p, search_multi, bacrun)

