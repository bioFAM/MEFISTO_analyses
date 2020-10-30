# Run this script from the terminal using: python downsample_evodevo.py Ndown seed species K

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
import os
import seaborn as sns
import matplotlib.pyplot as plt

def run_evodevo(nfactors = 5, Ndown = 3, warp = False, save = True,
              warping_ref = "Mouse", sample_seed = 4891,
              seed = 2020,
              species = ["Mouse", "Rabbit", "Rat", "Human", "Opossum"],
              views = ["Brain", "Cerebellum", "Heart", "Liver", "Testis"],
              model_groups = True, nm = None, tissue_as_sample = False):

    if tissue_as_sample:
        assert not warp, "Need to adapt warping reference if tissues are treated as groups"

    # specify data directory of normalized gene expression data
    if species == ["Mouse", "Rabbit", "Rat"] and not warp:
        nmtmp = "MRRab"
        datadir = "data/input_data/MRRab_matched/"
    elif warp:
        nmtmp = "warping"
        datadir = "data/input_data/all_unmatched/"
    else:
        print("Matched inputs are only provided for [Mouse, Rabbit, Rat]")
        sys.exit()


    # set filenames for output
    if nm is not None:
        nm = nm
    else:
        nm = nmtmp

    # load data and covariate
    data = []
    times = []
    samples_names = []
    if tissue_as_sample:
        group_names = []
        data_view = []
        for m in views:
            for g in species:
                df = pd.read_csv(datadir + "view_" + m + "_group_" + g + ".csv", header=0, index_col=0)
                data_view.append(np.asarray(df).transpose())
                times.append(np.asarray(
                    pd.read_csv(datadir + "times_group_" + g + ".csv", header=0, index_col=0)).transpose())
                samples_names.append(df.columns)
                group_names.append(m +"-"+g)
        data = [data_view]
        features_names = [df.index]
    else:
        for m in views:
            data_view = []
            for g in species:
                data_view.append(np.asarray(pd.read_csv(datadir + "view_" + m + "_group_" + g + ".csv", header = 0, index_col = 0)).transpose())
                if m == "Brain": # only needed once
                    times.append(np.asarray(pd.read_csv(datadir + "times_group_" + g + ".csv", header=0, index_col=0)).transpose())
            data.append(data_view)

    # convert warping ref to numeric
    warping_ref = np.where([species[i] == warping_ref for i in range(len(species))])[0][0]

    # mask values at random
    if Ndown > 0:
        np.random.seed(sample_seed)
        if tissue_as_sample:
            for i in range(len(data[0])):
                    Ng = data[0][i].shape[0]
                    masked_samples = np.random.choice(Ng, Ndown, replace=False)
                    data[0][i][masked_samples,:] = np.nan
        else:
            for m in range(len(views)):
                for g in range(len(species)):
                    Ng = data[m][g].shape[0]
                    masked_samples = np.random.choice(Ng, Ndown, replace=False)
                    data[m][g][masked_samples,:] = np.nan

    # check dimension and name views and groups
    if tissue_as_sample:
        assert len(data) == 1, "problem in loading data, wrong number of views"
        assert len(data[0]) == len(species) * len(views), "problem in loading data, wrong number of groups"
        view_names = ["mRNA"]
    else:
        assert len(data) == len(views), "problem in loading data, wrong number of views"
        assert len(data[0]) == len(species), "problem in loading data, wrong number of groups"
        view_names = views
        group_names = species

    # prepare MOFA model with time as covariate
    ent = entry_point()
    ent.set_data_options()
    ent.set_data_matrix(data, groups_names = group_names, views_names = view_names)
    ent.set_model_options(factors=nfactors)
    ent.set_train_options(seed=seed, convergence_mode = "medium")
    ent.set_covariates(times, covariates_names = "time")
    ent.set_smooth_options(warping=warp, warping_ref = warping_ref, model_groups = model_groups)

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
    if save:
        if Ndown == 0:
            if model_groups:
                outfile = "out/evodevo_groups_%s-seed_%s.hdf5" % (nm, seed)
            else:
                outfile = "out/evodevo_%s-seed_%s.hdf5" % (nm, seed)

            # interpolate for missing time points
            ent.predict_factor(new_covariates=ent.model.nodes["Sigma"].covariates)

        else:
            if model_groups:
                outfile = "out/evodevo_groups_%s-N%s-sample_seed_%s.hdf5" % (nm, Ndown, sample_seed)
            else:
                outfile = "out/evodevo_%s-N%s-sample_seed_%s.hdf5" % (nm, Ndown, sample_seed)

        ent.save(outfile)

    # write output to csv
    results = {'time': total, 'mem_usage': peak, 'n_down': Ndown, 'sample_seed' : sample_seed, 'seed' : seed}
    df = pd.DataFrame.from_dict(data=results, orient='index').T
    if model_groups:
        stats_file = 'out/evodevo_groups_%s_stats.csv' % nm
    else:
        stats_file = 'out/evodevo_%s_stats.csv' % nm
    if os.path.exists(stats_file):
        df.to_csv(stats_file, mode='a', header=False)
    else:
        df.to_csv(stats_file, header=True)


if __name__ == '__main__':

    Ndown =  sys.argv[1]
    sample_seed = int(sys.argv[2])
    K = sys.argv[4]

    if sys.argv[3] == "all":
        nm = 'all_k%s' % K
        run_evodevo(nfactors = K, Ndown = int(Ndown), sample_seed = sample_seed, warp = True,
         species = ["Mouse", "Rabbit", "Rat", "Human", "Opossum"],
         views = ["Brain", "Cerebellum", "Heart", "Liver", "Testis"], nm = nm)
    else:
        "Group specification not understood."