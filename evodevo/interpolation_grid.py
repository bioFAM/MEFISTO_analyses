# Script to run MOFA and MEFISTO for a comparison of interpolation results on evodevo data

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

# missing = view/sample combinatinos missing (full sample missingness in a view)
def run_grid(nfactors = 5, Tmissing = 5, GP_factors = True, note = "", masking_seed = 1234,
             Nviews = 1, seed = 2020, model_groups = True):

    nm = "all"
    species = ["Mouse", "Rabbit", "Rat", "Human", "Opossum"]
    views = ["Brain", "Cerebellum", "Heart", "Liver", "Testis"]

    # specify data directory of normalized gene expression data
    datadir = "data/input_data/all_unmatched/"

    # set and check number of views to mask a time point in
    if Nviews == "all":
        Nviews = len(views)

    if Nviews > len(views):
        print("Nviews is larger than available number of views, setting to all views.")
        Nviews = min(Nviews, len(views))

    # load data
    data = []
    times = []
    for m in views:
        data_view = []
        for g in species:
            data_view.append(np.asarray(pd.read_csv(datadir + "view_" + m + "_group_" + g + ".csv", header = 0, index_col = 0)).transpose())
            if m == "Brain": # only needed once
                times.append(np.asarray(pd.read_csv(datadir + "times_group_" + g + ".csv", header=0, index_col=0)).transpose())
        data.append(data_view)

    # check dimension
    assert len(data) == len(views), "problem in loading data, wrong number of groups"
    assert len(data[0]) == len(species), "problem in loading data, wrong number of views"

    # keep unmasked data in data_full and center as done in MOFA (per group)
    data_full = copy.deepcopy(data)
    for g in range(len(species)):
        for m in range(len(views)):
            data_full[m][g]-= np.nanmean(data_full[m][g], axis=0)

    # mask values (draw timepoint - species combinations and views at random)
    times_spec = np.vstack(
        [np.repeat(range(len(species)), [len(times[g]) for g in range(len(species))]),
        np.concatenate([np.arange(len(times[g])) for g in range(len(species))])]
    )
    np.random.seed(masking_seed)
    times_spec2mask = np.random.choice(range(times_spec.shape[1]), Tmissing, replace=False)
    for ts in times_spec2mask:
        g2mask = times_spec[0, ts]
        t2mask = times_spec[1, ts]
        views2mask = np.random.choice(range(len(views)), Nviews, replace=False)
        for m in views2mask:
            data[m][g2mask][t2mask,:] = np.nan

    # warping reference as in full training
    warping_ref = "Mouse"
    warping_ref = np.where([species[i] == warping_ref for i in range(len(species))])[0][0]


    # prepare MOFA model with/without covarites (determined by GP_factors)
    ent = entry_point()
    ent.set_data_options()
    ent.set_data_matrix(data, groups_names = species, views_names=views)
    ent.set_model_options(factors=nfactors)
    ent.set_train_options(seed=seed)
    if GP_factors:
        ent.set_covariates(times, covariates_names = "time")
        ent.set_smooth_options(warping=True, warping_ref = warping_ref, model_groups = model_groups)

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
    mse = 0
    n_missing = 0
    for m in range(len(views)):
        mse_m = np.nansum(((ent.imputed_data["mean"][m][ent.model.nodes['Y'].getNodes()[m].getMask()] - np.vstack(data_full[m])[ent.model.nodes['Y'].getNodes()[m].getMask()])**2)) # samples with NAs in full data are removed from sum, only consider masked ones
        mse = mse + mse_m
        n_missing = n_missing + ent.model.nodes['Y'].getNodes()[m].getMask().sum() - sum([np.isnan(data_full[m][g]).sum() for g in range(len(species))]) # count all missing values where data was available
    mse = mse / n_missing

    # write output to csv
    results = {'mse' : mse, 'time':  total, 'GP_factors' : GP_factors, 'n_missing': n_missing,
               'Tmissing' : Tmissing, 'masking_seed': masking_seed, 'date' : date.today(), 'note' : note, 'mem_usage': peak,
               'n_factors' : Zlearnt.shape[1], 'Nviews' : Nviews}
    df = pd.DataFrame.from_dict(data=results,orient='index').T

    if os.path.exists('out/interpol_results_evodevo.csv'):
        df.to_csv('out/interpol_results_evodevo.csv', mode='a', header=False)
    else:
        df.to_csv('out/interpol_results_evodevo.csv', header=True)


if __name__ == '__main__':

    seed = int(sys.argv[1])
    Tmissing = int(sys.argv[2])
    Nviews = int(sys.argv[3])
    # 
    run_grid(Tmissing = int(Tmissing), GP_factors = True, masking_seed = seed, Nviews = Nviews)
    run_grid(Tmissing = int(Tmissing), GP_factors = False, masking_seed = seed, Nviews = Nviews)

