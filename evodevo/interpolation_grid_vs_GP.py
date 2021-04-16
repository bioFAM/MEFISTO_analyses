# Script to run MOFA and MEFISTO for a comparison of interpolation results on evodevo data

import time
import numpy as np
import pandas as pd
import os
import sys
import tracemalloc
from datetime import date
import copy
sys.path.append("..")
from utils.interpolation import fit_GP, fit_MOFA, calc_mse

# missing = view/sample combinatinos missing (full sample missingness in a view)
def run_grid(nfactors = 5, Tmissing = 5, GP_factors = True, note = "", masking_seed = 1234,
             Nviews = 1, seed = 2020, model_groups = True, method = "FA", max_iter = 1000,
             species = ["Mouse", "Rabbit", "Rat", "Human", "Opossum"],
             views = ["Brain", "Cerebellum", "Heart", "Liver", "Testis"], nm = "all", warping = True,
             frac_features = 1, n_genes = 1000):


    M = len(views)
    G = len(species)

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
            dd_m_g = np.asarray(pd.read_csv(datadir + "view_" + m + "_group_" + g + ".csv", header = 0, index_col = 0)).transpose()
            data_view.append(dd_m_g)
            if m == views[0]: # only needed once
                times.append(np.asarray(pd.read_csv(datadir + "times_group_" + g + ".csv", header=0, index_col=0)).transpose())
        data.append(data_view)

    if n_genes != "all":
        np.random.seed(masking_seed + 2020)
        genes2keep = np.random.choice(range(data[0][0].shape[1]), n_genes, replace=False) # keep a subset of genes in all species and organs
        for m in range(M):
            for g in range(G):
                data[m][g] = data[m][g][:, genes2keep]

    # check dimension
    assert len(data) == M, "problem in loading data, wrong number of groups"
    assert len(data[0]) == G, "problem in loading data, wrong number of views"

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
    if frac_features < 1:
        D = data[0][0].shape[1]
        genes2mask = np.random.choice(range(D), int(frac_features * D), replace=False)
    for ts in times_spec2mask:
        g2mask = times_spec[0, ts]
        t2mask = times_spec[1, ts]
        views2mask = np.random.choice(range(len(views)), Nviews, replace=False)
        for m in views2mask:
            if frac_features < 1:
                data[m][g2mask][t2mask,genes2mask] = np.nan
            else:
                data[m][g2mask][t2mask, :] = np.nan

    data_masked = copy.deepcopy(data)

    # warping reference as in full training
    warping_ref = "Mouse"
    warping_ref = np.where([species[i] == warping_ref for i in range(len(species))])[0][0]

    tracemalloc.start()
    t0 = time.time()

    if method != "univGPs":
        predictions,ent = fit_MOFA(data = data,
                                   times = times,
                                   nfactors = nfactors,
                                   seed = seed,
                                   GP_factors = GP_factors,
                                   warping = warping,
                                   warping_ref = warping_ref,
                                   model_groups = model_groups,
                                   iter = max_iter)
    else:
        predictions, model, likelihood = fit_GP(data, times, iter = max_iter)
    t1 = time.time()
    total_time = t1 - t0
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # evaluate interpolation MSE (on missing values that hava a ground truth)
    mse, mse_mean, n_missing = calc_mse(data_masked, data_full, predictions)

    # write output to csv
    results = {'mse' : mse, 'mse_mean': mse_mean, 'time':  total_time, 'GP_factors' : GP_factors, 'n_missing': n_missing,
               'Tmissing' : Tmissing, 'masking_seed': masking_seed, 'date' : date.today(),
               'note' : note, 'mem_usage': peak, 'Nviews' : Nviews, 'method' : method, 'n_genes' : n_genes,
               'frac_features' : frac_features}
    df = pd.DataFrame.from_dict(data=results,orient='index').T

    filenm = 'out/interpol_results_evodevo_%s_%s.csv' % (nm, method)
    if os.path.exists(filenm):
        df.to_csv(filenm, mode='a', header=False)
    else:
        df.to_csv(filenm, header=True)


if __name__ == '__main__':

    seed = int(sys.argv[1])
    Tmissing = int(sys.argv[2])
    method = sys.argv[3]
    frac_features = float(sys.argv[4])


    if method != "univGPs":
        run_grid(Tmissing = int(Tmissing), GP_factors = True, masking_seed = seed,
                 species = ["Mouse"], views = ["Brain"], nm = "single_mod", warping=False, frac_features = frac_features)
        run_grid(Tmissing = int(Tmissing), GP_factors = False, masking_seed = seed,
                 species = ["Mouse"], views = ["Brain"], nm = "single_mod", warping=False, frac_features = frac_features)
    else:
        run_grid(Tmissing = int(Tmissing), GP_factors = True, masking_seed = seed, method = "univGPs",
                 species = ["Mouse"], views = ["Brain"], nm = "single_mod", warping=False, frac_features = frac_features)


