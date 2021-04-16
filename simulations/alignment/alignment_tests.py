# Script to run MOFA+ and MEFISTO for a comparison on simulated data
# across different # samples, # groups, noise levels and # missing values 

import time
import copy
import scipy.stats as ss
import numpy as np
import pandas as pd
import os
import sys
import tracemalloc
from datetime import date
from mofapy2.run.entry_point import entry_point
from mofapy2.simulate import simulate_mofa as simmofa
import matplotlib.pyplot as plt


# missing = view/sample combinations missing (full sample missingness)
def test_alignment(n_smooth_factors = 1, n_nonsmoothfactors = 3, nfactors_sim_extra = 0,  G = 3, N = 20, D = 500, noise_level = 1,
             seed = 1234567, M = 4, max_iter = 1000, verbose = False, save = False, missing = 0.1,
                   missing_all = 0.1, experiment = "test", weight_views = False):

    # simulate one smooth shared factor and K non-smooth, non-shared factors
    n_nonsmoothfactors = int(n_nonsmoothfactors)
    nfactors = n_smooth_factors + n_nonsmoothfactors
    scale_smooth = 1
    ls_smooth =  0.2
    sharedness = [True] * n_smooth_factors + [False] * n_nonsmoothfactors + [False] * nfactors_sim_extra
    lscales = [ls_smooth] * n_smooth_factors + [0] * n_nonsmoothfactors  + [0] * nfactors_sim_extra
    scales = [scale_smooth] * n_smooth_factors + [0] * n_nonsmoothfactors + [0] * nfactors_sim_extra
    groupsidx = np.repeat(range(G), N)

    # simulate data
    np.random.seed(seed)
    if experiment == "imbalance" or experiment == "imbalance_weighted":
        assert M == 2 and n_smooth_factors == 1 and n_nonsmoothfactors == 1
        inactive = 1000
        active = 1
        # in view 1 first (smooth) factor is active, 2 (non-smooth) inactive; reverse in view 2
        alpha = [np.array([active] +[inactive]  + [inactive]*nfactors_sim_extra), np.array([inactive, active]  + [active]*nfactors_sim_extra)]
    else:
        alpha = None

    if isinstance(D, int):
        D = [D] * M
    else:
        assert isinstance(D, list) and len(D) == M

    sim = simmofa.simulate_data(N=N, seed=seed, views=list(range(M)), D=D,
                                K=nfactors + nfactors_sim_extra, G=G, lscales=lscales, noise_level=noise_level,
                                scales = scales, shared = sharedness, alpha = alpha)

    # mask parts of the data
    data_full = copy.deepcopy(sim['data'])
    sim['data'] = simmofa.mask_samples(sim, perc = missing, perc_all_views = missing_all)

    # misalign covariates between groups
    sim['sample_cov'][1] = np.exp(sim['sample_cov'][1])
    sim['sample_cov'][2] = 0.4 * sim['sample_cov'][2] + 0.3
    

    # prepare model
    ent = entry_point()
    ent.set_data_options(scale_views=False)
    ent.set_data_matrix(sim['data'])
    ent.set_model_options(factors=nfactors)
    ent.set_train_options(seed=2020, convergence_mode="fast", iter=max_iter, verbose=verbose, weight_views = weight_views)
    ent.set_covariates(sim['sample_cov'])
    ent.set_smooth_options(warping=True, model_groups = True)

    # run and build the model
    tracemalloc.start()
    ent.build()
    t0 = time.time()
    ent.run()
    t1 = time.time()
    total = t1 - t0
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # get inferred hyperparameters
    scales_learnt = ent.model.train_stats['scales']
    lscales_learnt = ent.model.train_stats['length_scales']

    # get factors recovery error
    Zlearnt = ent.model.getExpectations()['Z']['E']

    # calculate factor recovery error    
    # if not right number of factors inferred set to nan
    if not Zlearnt.shape[1] == nfactors: 
        factor_r2 = np.nan
        scales_learnt = [np.nan] * nfactors
        lscales_learnt = [np.nan] * nfactors
        factor_idx = [np.nan] * nfactors
    else:
        Zsim = np.vstack(sim['Z'])
         # get idx of learnt factor corresponding to simulated factor by maximal correlation
        factor_idx = [np.argmax([abs(ss.pearsonr(Zsim[:,p], Zlearnt[:,pp])[0]) for pp in range(nfactors)]) for p in range(nfactors)]
        #check for duplicates - if true not all facotrs are captured on a unique factor
        if not len(factor_idx) == len(set(factor_idx)): 
            factor_r2 = np.nan
        else:
            # calculate correlation between inferred and simulated factors 
            factor_r2 = np.mean([np.max([abs(ss.pearsonr(Zsim[:, pp], Zlearnt[:, p])[0]) for pp in range(nfactors)]) for p in
         range(nfactors)]) ** 2

            scales_learnt = scales_learnt[factor_idx] # match to simulated factor
            lscales_learnt = lscales_learnt[factor_idx]
            if verbose:
                print(scales_learnt)

    # get imputation error
    ent.impute(mask_outliers = False)
    mse = 0
    n_missing = 0
    imp_r2 = 1
    if missing + missing_all > 0:
        for m in range(M):
            mse_m = ((ent.imputed_data["mean"][m][ent.model.nodes['Y'].getNodes()[m].getMask()] - np.vstack(data_full[m])[ent.model.nodes['Y'].getNodes()[m].getMask()])**2).sum()
            mse = mse + mse_m
            n_missing = n_missing + ent.model.nodes['Y'].getNodes()[m].getMask().sum()
        mse = mse / n_missing

        imp_r2 = np.mean([ss.pearsonr(np.vstack(data_full[m])[ent.model.nodes['Y'].getNodes()[m].getMask()].flatten(), ent.imputed_data["mean"][m][ent.model.nodes['Y'].getNodes()[m].getMask()].flatten())[0] ** 2 for m in range(M)]) 
    rec_r2 = np.mean([ss.pearsonr(np.vstack(data_full[m]).flatten(), ent.imputed_data["mean"][m].flatten())[0]  ** 2 for m in range(M)] )

    # get warping error
    sample_cov_transformed = ent.model.getNodes()['Sigma'].sample_cov_transformed

    # compare to untransformed group
    warp_mse = sum([sum((sample_cov_transformed[groupsidx == g] - sim['sample_cov'][0])**2) for g in range(G)]) / (N * G)
    nowarp_mse = sum([sum((sim['sample_cov'][g] - sim['sample_cov'][0])**2) for g in range(G)]) / (N * G)
    
    warp_mse = warp_mse[0]
    nowarp_mse = nowarp_mse[0]

    warp_score = warp_mse / nowarp_mse

    # get group covariance error:
    if len(factor_idx) == len(set(factor_idx)):
        Gmat_learnt = ent.model.nodes['Sigma'].getParameters()['Kg']

        # get sharedness error
        true_sharedness = [np.mean(np.abs(sim['Gmats'][k] - np.eye(G))[np.triu_indices(G, 1)]) for k in range(nfactors)]
        inferred_sharedness = [np.mean(np.abs(Gmat_learnt[factor_idx[k]] - np.eye(G))[np.triu_indices(G, 1)]) for k in
                               range(nfactors)]

    # if no group covariance was simulated set to nan     
    else:
        true_sharedness = [np.nan] * nfactors
        inferred_sharedness = [np.nan] * nfactors

    # write output to csv
    results = {'factor_r2' : factor_r2, 'time':  total, 'n_nonsmoothfactors' : n_nonsmoothfactors,
                'N' : N, 'G': G, 'D1' : D[1],  'D0' : D[0], 'weight_views' : weight_views,
               'D' : D, 'noise_level': noise_level, 'missing' : missing + missing_all, 'seed' : seed,
               'date' : date.today(), 'mem_usage': peak, 'nfactors_sim_extra' : nfactors_sim_extra,
               'n_factors' : nfactors,
               'warp_mse' : warp_mse,
               'warp_score' : warp_score,
               'nowarp_mse' : nowarp_mse,
               'mse' : mse,  'imp_r2' : imp_r2, 'rec_r2' : rec_r2}
    if verbose:
        print(results)

    df = pd.DataFrame.from_dict(data=results, orient='index').T

    # save summary statistics if not the model itself is saved
    outfile = "out/" + experiment + ".csv"
    if not save:
        if os.path.exists(outfile):
            df.to_csv(outfile, mode='a', header=False)
        else:
            df.to_csv(outfile, header=True)

    else:
        ent.save("out/" + experiment +".hdf5")


if __name__ == '__main__':

    # seed for simulations
    seed = int(sys.argv[1])
    experiment = sys.argv[2]

    if experiment == "nonsmooth":
        for K_nonsmooth in [15, 12, 9, 6, 3, 2, 1, 0]:
            test_alignment(n_nonsmoothfactors = K_nonsmooth, seed=seed, experiment = experiment)
        test_alignment(n_nonsmoothfactors=15, n_smooth_factors=0, seed=seed, experiment=experiment)
    elif experiment == "imbalance":
        for D_large in [100, 200, 500, 1000, 2000, 3000, 4000, 5000]:
            test_alignment(n_smooth_factors=1, n_nonsmoothfactors=1, seed=seed, experiment=experiment,
                           M = 2, D=[100, D_large])
    elif experiment == "imbalance_weighted":
        for D_large in [100, 200, 500, 1000, 2000, 3000, 4000, 5000]:
            test_alignment(n_smooth_factors=1, n_nonsmoothfactors=1, seed=seed, experiment=experiment,
                           M = 2, D=[100, D_large], weight_views = True)
    # elif experiment == "imbalance_extra_weighted":
    #     for D_large in [100, 200, 500, 1000, 2000, 3000, 4000, 5000]:
    #         print("Use weighting")
    #         test_alignment(n_smooth_factors=1, n_nonsmoothfactors=1, seed=seed, experiment=experiment,
    #                        M=2, D=[100, D_large], weight_views = True, nfactors_sim_extra = 2)
    #         test_alignment(n_smooth_factors=1, n_nonsmoothfactors=1, seed=seed, experiment=experiment,
    #                        M=2, D=[100, D_large], weight_views = True, nfactors_sim_extra = 1)
    # elif experiment == "imbalance_extra":
    #     for D_large in [100, 200, 500, 1000, 2000, 3000, 4000, 5000]:
    #         print("No weighting")
    #         test_alignment(n_smooth_factors=1, n_nonsmoothfactors=1, seed=seed, experiment=experiment,
    #                        M=2, D=[100, D_large], weight_views = False, nfactors_sim_extra = 2)
    #         test_alignment(n_smooth_factors=1, n_nonsmoothfactors=1, seed=seed, experiment=experiment,
    #                        M=2, D=[100, D_large], weight_views = False, nfactors_sim_extra = 1)
    else:
        print("Experiment not recognised.")