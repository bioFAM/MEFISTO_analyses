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
def run_grid(nfactors = 3, G = 5, N = 20, Dm = 500, noise_level = 1, missing = 0.1, missing_all = 0.1, seed = 1234567,
              method = "MEFISTO", note = "none", lscales = [0.2, 0.1, 0.0], scales = [1, 0.6, 0], M = 4, plot = False,
             max_iter = 1000, verbose = False, sparse_frac = 0.75, warp = False, save = False, group_differences = True,
             model_groups = True):

    nfactors = int(nfactors)
    assert len(lscales) == nfactors
    assert len(scales) == nfactors

    groupsidx = np.repeat(range(G), N)

    # simulate data
    np.random.seed(seed)
    if group_differences:
        if nfactors == 3:
            sharedness = np.random.choice([True, False], 2, replace=False).tolist() + [False] # one shared, one non-shared, one non-smooth
        else:
            sharedness = True # not in use
    else:
        sharedness = True
    sim = simmofa.simulate_data(N=N, seed=seed, views=["0", "1", "2", "3"], D=[Dm] * M,
                                K=nfactors, G=G, lscales=lscales, noise_level=noise_level,
                                scales = scales, shared = sharedness)

    # mask parts of the data
    data_full = copy.deepcopy(sim['data'])
    sim['data'] = simmofa.mask_samples(sim, perc = missing, perc_all_views = missing_all)

    # misalign covariates between groups
    if warp:
        assert G == 3, "Warping defined only for G=3"
        sim['sample_cov'][1] = np.exp(sim['sample_cov'][1])
        sim['sample_cov'][2] = 0.4 * sim['sample_cov'][2] + 0.3
    
    # optional plotting of simulated factors
    if plot:
        fig, axs = plt.subplots(1, nfactors)
        Zsim = sim['Z']
        for g in range(G):
            for i in range(nfactors):
                axs[i].scatter(sim['sample_cov'][g], Zsim[g][:, i])
                axs[i].set_title("simulated factors")

    # prepare model
    ent = entry_point()
    ent.set_data_options(scale_views=False)
    ent.set_data_matrix(sim['data'])
    ent.set_model_options(factors=nfactors)
    ent.set_train_options(seed=2020, convergence_mode="fast", iter=max_iter, verbose=verbose)

    # for time-aware multi-modal FA with GP model add covariates
    if not method == "MOFA2":
        ent.set_covariates(sim['sample_cov'])
        if method == "MEFISTO+align":
            ent.set_smooth_options(warping=True, model_groups = model_groups)
        elif method == "MEFISTO_sparse":
            ent.set_smooth_options(model_groups = model_groups, sparseGP = True, n_inducing= int((N * G) * sparse_frac))
        else:
            ent.set_smooth_options(model_groups = model_groups)


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
    if method != "MOFA2":
        scales_learnt = ent.model.train_stats['scales']
        lscales_learnt = ent.model.train_stats['length_scales']
    else:
        scales_learnt = np.array([np.nan] * nfactors)
        lscales_learnt = np.array([np.nan] * nfactors)

    # get factors recovery error
    Zlearnt = ent.model.getExpectations()['Z']['E']

    # calculate factor recovery error    
    # if not right number of factors inferred set to nan
    if not Zlearnt.shape[1] == nfactors: 
        factor_r2 = np.nan
        scales_learnt = [np.nan] * nfactors
        lscales_learnt = [np.nan] * nfactors
        post_var = [np.nan] * nfactors
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
            if method != "MOFA2":
                scales_learnt = scales_learnt[factor_idx] # match to simulated factor
                lscales_learnt = lscales_learnt[factor_idx]
                if verbose:
                    print(scales_learnt)
        
        # get posterior variance
        post_var = ent.model.getExpectations()['Z']['E2'] - (ent.model.getExpectations()['Z']['E']) ** 2
        post_var = post_var.mean(axis = 0)

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
    if method == "MEFISTO+align":
        sample_cov_transformed = ent.model.getNodes()['Sigma'].sample_cov_transformed
        # compare to untransformed group
        warp_mse = sum([sum((sample_cov_transformed[groupsidx == g] - sim['sample_cov'][0])**2) for g in range(G)]) / (N * G)

    else: # no transformation made
        warp_mse = sum([sum((sim['sample_cov'][g] - sim['sample_cov'][0])**2) for g in range(G)]) / (N * G)
   
    warp_mse = warp_mse[0]

    # get group covariance error:
    if group_differences and len(factor_idx) == len(set(factor_idx)):
        if "Sigma" in ent.model.nodes.keys() and model_groups:
            Gmat_learnt = ent.model.nodes['Sigma'].getParameters()['Kg']
        # MEFISTO without model_groups assumes all groups to be connected
        elif "Sigma" in ent.model.nodes.keys():
            Gmat_learnt = [np.ones([G,G])] * nfactors
        # MOFA2 assumes all groups to be unconnected
        else:
            Gmat_learnt = [np.eye(G)] * nfactors

        # get sharedness error
        true_sharedness = [np.mean(np.abs(sim['Gmats'][k] - np.eye(G))[np.triu_indices(G, 1)]) for k in range(nfactors)]
        inferred_sharedness = [np.mean(np.abs(Gmat_learnt[factor_idx[k]] - np.eye(G))[np.triu_indices(G, 1)]) for k in
                               range(nfactors)]

    # if no group covariance was simulated set to nan     
    else:
        true_sharedness = [np.nan] * nfactors
        inferred_sharedness = [np.nan] * nfactors

    # write output to csv
    results = {'factor_r2' : factor_r2, 'time':  total, 'method' : method, 'model_groups' : model_groups,
                'group_differences': group_differences, 'N' : N, 'G': G,
               'Dm' : Dm, 'noise_level': noise_level, 'missing' : missing + missing_all, 'seed' : seed,
               'date' : date.today(), 'note' : note, 'mem_usage': peak, 'lscales' : lscales,
               'scales' : scales, 'sparse_frac' : sparse_frac,
               'n_factors' : nfactors,
               'warp_mse' : warp_mse,
               'n_factors_learnt' : Zlearnt.shape[1],
               'scales_learnt' : scales_learnt,
               'lscales_learnt' : lscales_learnt,
               'true_sharedness' : np.array(true_sharedness),
               'inferred_sharedness' : np.array(inferred_sharedness),
               'post_var' : post_var, 'mse' : mse,  'imp_r2' : imp_r2, 'rec_r2' : rec_r2}
    if verbose:
        print(results)

    df = pd.DataFrame.from_dict(data=results, orient='index').T
    # expand multi-factor columns
    for nm in ['scales', 'lscales', 'scales_learnt', 'lscales_learnt', 'true_sharedness', 'inferred_sharedness', 'post_var']:
        dfsplit = df[nm].apply(pd.Series)
        dfsplit = dfsplit.rename(columns=lambda x: nm + "_" + str(x))
        df = pd.concat([df, dfsplit], axis=1)
        df = df.drop(columns = [nm])

    # optional plotting of inferred factors
    if plot:
        Zlearnt = ent.model.getExpectations()['Z']['E']
        fig, axs = plt.subplots(1, nfactors)
        for g in range(G):
            for i in range(nfactors):
                axs[i].scatter(sim['sample_cov'][g], Zlearnt[groupsidx == g, i])
                axs[i].set_title("inferred factors")

    # save summary statistics if not the model itself is saved
    if not save:
        if os.path.exists('out/simulation_results.csv'):
            df.to_csv('out/simulation_results.csv', mode='a', header=False)
        else:
            df.to_csv('out/simulation_results.csv', header=True)

    else:
        ent.save("out/grid_model.hdf5")


if __name__ == '__main__':
    # seed for simulations
    seed = int(sys.argv[2])
    method = sys.argv[3]  

    # run different experiments depending on first argument
    # 1 - sparse
    if sys.argv[1] == "sparse":
        for N in np.linspace(500,2500,5): 
            run_grid(N=int(N), G = 1,  method=method, note="sparse", seed=seed)

    # 2 - align
    elif sys.argv[1] == "align":
        for N in np.linspace(10, 100, 10):
            run_grid(N=int(N), G = 3, method=method, note="align", seed=seed, warp = True)

    # 2 - illustrative example model for alignemt (with one non-smooth and two shared factors)
    elif sys.argv[1] == "align_example":
        run_grid(N = 20, G = 3, method = "MEFISTO+align", note="align_example", seed=seed, warp = True,
         save = True, group_differences = False)

    # 3- base grid: varying noise, missing, N and G
    else:
        # number of samples per group
        if sys.argv[1] == "N":
            for N in np.linspace(10, 100, 10):
                run_grid(N =int(N), method = method, note = "vary_N", seed = seed)
        # number of groups
        elif sys.argv[1] == "G":
            for G in np.linspace(1, 9, 5):
                run_grid(G =int(G), method = method, note = "vary_G", seed = seed)
        # noise level
        elif sys.argv[1] == "noise":
            for noise_level in [0.01, 0.1, 0.5, 1, 2, 3]:
                run_grid(noise_level=noise_level, method=method, note = "vary_noise", seed = seed)
        # missing values
        elif sys.argv[1] == "missing":
            for missing in np.linspace(0,0.45,10): # note: the total level is 2 x this
                run_grid(missing=missing, missing_all = missing, method=method, note = "vary_missing", seed = seed)
        else:
            print("Experiment name not recognised.")
