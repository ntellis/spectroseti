import matplotlib
from pathos.multiprocessing import ProcessingPool as Pool

matplotlib.use('TkAgg')

import spectroseti.runner as runner

observations = [['awx', 222]]  # ,[],['awx',223],['awx',224]]

LS = runner.LaserSearch()

#LS.search_multiple(observations, output_pngs=1, number_mads=8)

p = Pool(6)
LS.search_multiple(observations, output_pngs=1,number_mads=5)


