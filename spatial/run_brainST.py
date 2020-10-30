# script to run MEFISTO on spatial transcriptomics data
# call this script using python run_brainST.py n_inducing n_factors

from mofapy2.run.entry_point import entry_point

import pandas as pd
import scipy as s
import numpy as np
import io
import h5py
import scipy.stats as stats
import matplotlib.pyplot as plt
import math
import scipy.spatial as SS
import pandas as pd
import os
import sys
import tracemalloc
import time


# load data
dd = pd.read_csv("data/brain_view_1.csv", index_col=0)
scov = pd.read_csv("data/brain_sample_cov.csv", index_col=0)
dd = np.asarray(dd).transpose()
scov = np.asarray(scov).transpose()

###########################
## Initialise MEFISTO model ##
###########################

# initialise the entry point
ent = entry_point()

# Set data and options
n_factors = sys.argv[2]
ent.set_data_options(scale_views=False)
ent.set_data_matrix([[dd]])
ent.set_model_options(factors = int(n_factors))
ent.set_train_options(seed = 2020)

# Set covariate and smooth options
ent.set_covariates([scov])
n_inducing = sys.argv[1]
if n_inducing != "full":
	ent.set_smooth_options(sparseGP = True, n_inducing = int(n_inducing))
else:
	ent.set_smooth_options(sparseGP = False)


# Build and run the model
tracemalloc.start()
ent.build()
t0 = time.time()
ent.run()
t1 = time.time()
total = t1 - t0
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()


# save model
outfile= "out/brain_N%s_K%s.hdf5" % (n_inducing, n_factors)
ent.save(outfile)

# write output to csv
results = {'time': total, 'mem_usage': peak, 'n_inducing' : n_inducing, 'n_factors' : n_factors}
df = pd.DataFrame.from_dict(data=results, orient='index').T
if os.path.exists('out/brain_stats.csv'):
    df.to_csv('out/brain_stats.csv', mode='a', header=False)
else:
    df.to_csv('out/brain_stats.csv', header=True)