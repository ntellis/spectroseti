import numpy as np

no_cosmic = np.load('scripts/no_cosmic_examples.npy')
cosmic = np.load('scripts/cosmic_examples.npy')

all_cosmic_ud_ratio = np.abs(cosmic[1:,:,1:]/cosmic[:-1,:,1:] - 1)
all_cosmic_lr_ratio = np.abs(cosmic[:,1:,1:]/cosmic[:,:-1,1:] - 1)
all_lr_ratio = np.abs(no_cosmic[:,1:,1:]/no_cosmic[:,:-1,1:] - 1)
all_ud_ratio = np.abs(no_cosmic[1:,:,1:]/no_cosmic[:-1,:,1:] - 1)



sb.distplot(np.max(all_cosmic_lr_ratio.reshape((all_cosmic_lr_ratio.shape[0]*all_cosmic_lr_ratio.shape[1],-1)), axis=0), bins=30)



