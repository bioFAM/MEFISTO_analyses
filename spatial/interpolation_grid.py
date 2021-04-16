# Script to run MOFA and MEFISTO for a comparison of interpolation results on spatial data

import time
import itertools
import sklearn
import scipy.stats as ss
import numpy as np
import pandas as pd
import os
import sys
import tracemalloc
from datetime import date
from mofapy2.run.entry_point import entry_point
import copy
import mofapy2
if not '/icgc/dkfzlsdf/' in os.getcwd():
    sys.path.append('spatial/')

def run_grid(nfactors = 4, spots_missing = 100, GP_factors = True, note = "", masking_seed = 1234,
             seed = 2020, frac_features = 1, n_inducing = 1000):
    
    # load data
    dd = pd.read_csv("data/brain_view_1.csv", index_col=0)
    scov = pd.read_csv("data/brain_sample_cov.csv", index_col=0)
    dd = np.asarray(dd).transpose()
    scov = np.asarray(scov).transpose()
    n_spots = dd.shape[0]
    n_genes = dd.shape[1]

    # keep unmasked data in data_full and center
    data_full = copy.deepcopy(dd)
    data_full-= np.nanmean(data_full, axis=0)

    # mask values
    np.random.seed(masking_seed)
    spots2mask = np.random.choice(range(n_spots), spots_missing, replace=False)
    if frac_features < 1:
        for j in spots2mask:
            genes2mask = np.random.choice(range(n_genes), int(frac_features * n_genes), replace=False)
            dd[j, genes2mask] = np.nan
    else:
        dd[spots2mask, :] = np.nan

    # initialise the entry point
    ent = entry_point()

    # Set data and options
    ent.set_data_options(scale_views=False)
    ent.set_data_matrix([[dd]])
    ent.set_model_options(factors=int(nfactors))
    ent.set_train_options(seed=seed)

    # Set covariate and smooth options
    if GP_factors:
        ent.set_covariates([scov])
        if n_inducing < n_spots:
            ent.set_smooth_options(sparseGP=True, frac_inducing= n_inducing/n_spots)
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

    # make predictions
    Zlearnt = ent.model.nodes['Z'].getExpectation()
    ent.impute(mask_outliers = False)
    masked = ent.model.nodes['Y'].getNodes()[0].getMask()
    mse = np.nansum(((ent.imputed_data["mean"][0][masked] - data_full[masked])**2)) / masked.sum()
    print(mse)
    # write output to csv
    results = {'mse' : mse, 'time':  total, 'GP_factors' : GP_factors, 'spots_missing': spots_missing,
               'masking_seed': masking_seed, 'date' : date.today(), 'frac_features' : frac_features, 'mem_usage': peak,
               'n_factors' : Zlearnt.shape[1], 'n_inducing' : n_inducing}
    df = pd.DataFrame.from_dict(data=results,orient='index').T

    if os.path.exists('out/interpol_results_spatial.csv'):
        df.to_csv('out/interpol_results_spatial.csv', mode='a', header=False)
    else:
        df.to_csv('out/interpol_results_spatial.csv', header=True)


if __name__ == '__main__':

    seed = int(sys.argv[1])
    spots_missing = int(sys.argv[2])
    n_inducing = int(sys.argv[3])

    run_grid(spots_missing = spots_missing, GP_factors = True, masking_seed = seed, n_inducing = n_inducing)
    run_grid(spots_missing = spots_missing, GP_factors = False, masking_seed = seed, n_inducing = n_inducing)

