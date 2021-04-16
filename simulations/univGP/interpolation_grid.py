# Script to run MOFA+ and MEFISTO and univariate GPs for a interpolation comparison on simulated data
# across different # samples, # groups, noise levels and # missing values

import time
import numpy as np
import pandas as pd
import tracemalloc
import copy
import os
import sys
from datetime import date
from mofapy2.simulate import simulate_mofa as simmofa
sys.path.append("../../../")
from utils.interpolation import fit_MOFA, fit_GP, calc_mse, save_predictions


def run_grid(nfactors = 2, G = 1, N = 20, Dm = 100, M = 4,
             noise_level = 1, missing = 0.1,
             missing_all = 0.1, seed = 1234567,
             method = "MEFISTO", note = "none",
             lscales = [0.2, 0.1], scales = [1, 0.6], n_factors_sim = None,
             save = False, max_iter = 1000):

    nfactors = int(nfactors)
    if n_factors_sim is None:
        n_factors_sim = nfactors

    assert len(lscales) == n_factors_sim
    assert len(scales) == n_factors_sim

    groupsidx = np.repeat(range(G), N)

    # simulate data
    np.random.seed(seed)
    views_nms = [str(m) for m in range(M)]
    sim = simmofa.simulate_data(N=N, seed=seed, views=views_nms, D=[Dm] * M,
                                K=n_factors_sim, G=G, lscales=lscales, noise_level=noise_level,
                                scales = scales, shared = True)

    # keep unmasked data in data_full and center as done in MOFA (per group)
    data_full = copy.deepcopy(sim['data'])
    for g in range(G):
        for m in range(M):
            data_full[m][g]-= np.nanmean(data_full[m][g], axis=0)


    # mask parts of the data
    sim['data'] = simmofa.mask_samples(sim, perc = missing, perc_all_views = missing_all)
    data_masked = copy.deepcopy(sim['data'])

    # fit models
    tracemalloc.start()
    t0 = time.time()

    if method != "univGPs":
        if method == "MEFISTO":
            GP_factors = True
        else:
            GP_factors = False

        predictions, ent = fit_MOFA(data = sim['data'],
                                   times = sim['sample_cov'],
                                   nfactors = nfactors,
                                   seed = 2020,
                                   GP_factors = GP_factors,
                                   warping = False,
                                   warping_ref = 0,
                                   model_groups = True,
                                   iter = max_iter)

    else:
        predictions, model, likelihood = fit_GP(data = sim['data'],
                                                times = sim['sample_cov'],
                                                iter = max_iter)

    t1 = time.time()
    total_time = t1 - t0
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # evaluate interpolation MSE
    mse, mse_mean, n_missing = calc_mse(data_masked = data_masked,
                                        data_full = data_full,
                                        predictions = predictions)

    # save results
    outdir = 'out/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # save summary statistics if not the model itself is saved
    if not save:
        results = {'time': total_time, 'method': method,
                   'N': N, 'G': G, 'M': M, 'Dm': Dm,
                   'noise_level': noise_level, 'missing': missing, 'missing_all': missing_all,
                   'seed': seed, 'date': date.today(), 'note': note,
                   'mem_usage': peak, 'scales' : scales, 'lscales' :lscales,
                   'n_factors': nfactors, 'n_factors_sim' : n_factors_sim,
                   'mse': mse, 'mse_mean': mse_mean,
                   'n_missing': n_missing}

        df = pd.DataFrame.from_dict(data=results, orient='index').T
        for nm in ['scales', 'lscales']:  # expand multi-factor columns
            dfsplit = df[nm].apply(pd.Series)
            dfsplit = dfsplit.rename(columns=lambda x: nm + "_" + str(x))
            df = pd.concat([df, dfsplit], axis=1)
            df = df.drop(columns=[nm])

        filenm = 'interpolation_results_%s.csv' % method
        if os.path.exists(outdir + filenm):
            df.to_csv(outdir + filenm, mode='a', header=False)
        else:
            df.to_csv(outdir + filenm, header=True)

    # save model + predictions
    else:
        if method != "univGPs":
            ent.save(outdir + "grid_model.hdf5")

        save_predictions(predictions = predictions,
                         data_masked = data_masked,
                         times = sim['sample_cov'],
                         method = method,
                         outdir= outdir)

        save_predictions(predictions = [np.vstack(data_full[m]) for m in range(M)],
                         data_masked = data_masked,
                         times = sim['sample_cov'],
                         method = "ground_truth",
                         outdir = outdir)



if __name__ == '__main__':
    if sys.argv[1] == "--mode=client":
        for method in ["MEFISTO", "univGPs", "MOFA2"]:
            run_grid(method=method, note="baseline", seed=93243,
                     scales=[0], missing_all=0.3,
                     lscales=[0.2],
                     missing=0, nfactors=1, save = True)
    else:
        # seed for simulations
        print(sys.argv)
        seed = int(sys.argv[2])
        method = sys.argv[3]

        # #test
        if sys.argv[1] == "test":
            pass
        #     print("Testing mode...")
        #     run_grid(scales = [0, 0], method = method,
        #              note = "test", seed = seed, save = True)

        # baseline - should be the same for all methods
        # elif sys.argv[1] == "baseline":
        #     for missing in np.linspace(0.1, 0.9, 9):
        #         run_grid(method = method, note = "baseline", seed = seed,
        #                  scales = [0], missing_all=missing,
        #                  lscales = [0.2],
        #                  missing = 0, nfactors = 1)
        # number of timepoints
        elif sys.argv[1] == "N":
            for N in np.linspace(10, 50, 5):
                run_grid(N =int(N), method = method, note = "vary_N", seed = seed)
        # number of features
        elif sys.argv[1] == "D":
            for D in [25, 50, 75, 100, 125, 150, 175, 200]:
                run_grid(Dm =int(D), method = method, note = "vary_Dm", seed = seed)
        # missing time points
        elif sys.argv[1] == "missing":
            for missing in np.linspace(0.05,0.45,9):
                run_grid(missing=missing, missing_all = missing, method=method, note = "vary_missing", seed = seed)
        # smoothness
        elif sys.argv[1] == "smoothness":
            for scale2 in np.linspace(0,1,11):
                run_grid(scales = [scale2, scale2], method=method, note = "vary_smoothness", seed = seed) # previouse run 0.6 for one
        elif sys.argv[1] == "n_factors":
            for n_factors_sim in np.arange(2,11):
                run_grid(n_factors_sim = n_factors_sim, scales=[0.8]*n_factors_sim,
                         lscales=[0.2]*n_factors_sim, method=method, note="vary_n_factors_sim", seed=seed)
        else:
            print("Experiment name not recognised.")
