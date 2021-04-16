import pickle
import os
import numpy as np
import pandas as pd
import re
import sys


def load_data(datadir, metadir):
    # load data and covariate
    samples_meta = pd.read_csv(metadir + "samples_metadata.csv")
    child_ids = list(set(samples_meta.child_id))
    data = []
    times = []
    sample_names = []
    view_names = ["microbiome"]
    for m in view_names:
        data_view = []
        for g in child_ids:
            dd = pd.read_csv(datadir + "view_" + m + "_group_" + g + ".csv", header = 0, index_col = 0)
            samp_names = dd.columns.tolist()
            if g == child_ids[0]:
                feature_names = dd.index.tolist()
            data_view.append(np.asarray(dd).transpose())
            sample_names.append(samp_names)
            times.append(np.asarray(pd.read_csv(datadir + "times_group_" + g + ".csv", header=0, index_col=0)).transpose())
        data.append(data_view)
    group_names = child_ids
    feature_names = [feature_names]

    return data, times, group_names, sample_names, feature_names, view_names

def load_data_views(datadir, metadir):
    # load data and covariate
    child_ids = pd.read_csv(datadir + "groups.csv").x.tolist()
    view_names = pd.read_csv(datadir + "views.csv").x.tolist()
    data = []
    times = []
    sample_names = []
    feature_names = []
    for m in view_names:
        data_view = []
        for g in child_ids:
            dd = pd.read_csv(datadir + "view_" + m + "_group_" + g + ".csv", header = 0, index_col = 0)
            if m == view_names[0]:
                samp_names = dd.columns.tolist()
                sample_names.append(samp_names)
                times.append(np.asarray(pd.read_csv(datadir + "times_group_" + g + ".csv", header=0, index_col=0)).transpose())
            if g == child_ids[0]:
                feature_names.append(dd.index.tolist())
            data_view.append(np.asarray(dd).transpose())
        data.append(data_view)
    group_names = child_ids

    return data, times, group_names, sample_names, feature_names, view_names


def train_MOFA(input_data, times, group_names, feature_names, sample_names, view_names, outfile, use_GP = True, model_groups = True, center_groups = False):

    from mofapy2.run.entry_point import entry_point

    # prepare MEFISTO model
    ent = entry_point()
    ent.set_data_options(center_groups=center_groups)
    ent.set_data_matrix(input_data, groups_names=group_names,
                        features_names=feature_names, samples_names=sample_names,
                        views_names=view_names)

    ent.set_model_options(factors=2)
    ent.set_train_options(seed=2020)

    if use_GP:
        ent.set_covariates(times, covariates_names="month")
        ent.set_smooth_options(model_groups=model_groups, warping=False, warping_ref=0, n_grid=10, opt_freq=50, start_opt = 50) # opt_freq added for RCLR, set

    # Build and run the model
    ent.build()
    ent.run()

    # interpolate
    if use_GP:
        ent.predict_factor(new_covariates=ent.model.nodes["Sigma"].covariates)
    
    ent.save(outfile )


def train_NMF(input_data, feature_names, sample_names, outfile):
    from sklearn.decomposition import NMF

    X = []
    for m in range(len(data)):
        X.append(np.vstack(data[m]))  # concatenate samples across groups
    X = np.hstack(X) # concatenate features across views

    keep_sample = (~np.isnan(X)).sum(axis=1) > 0
    sample_names = np.concatenate(sample_names)[keep_sample]
    X = X[keep_sample, :]

    model = NMF(n_components=2, init='random', random_state=0)
    W = model.fit_transform(X)
    H = model.components_

    # save
    pd.DataFrame(W, index=sample_names).to_csv(outfile +"_NMF_W.csv")
    pd.DataFrame(H.transpose(), index=np.concatenate(feature_names)).to_csv(outfile +"_NMF_H.csv")


def train_ZIFA(input_data, feature_names, sample_names, outfile, use_block = False):
    from ZIFA import ZIFA

    X = []
    for m in range(len(data)):
        X.append(np.vstack(data[m]))  # concatenate samples across groups
    X = np.hstack(X) # concatenate features across views

    keep_sample = (~np.isnan(X)).sum(axis=1) > 0
    sample_names = np.concatenate(sample_names)[keep_sample]
    X = X[keep_sample, :]

    if not use_block:
        Z, model_params = ZIFA.fitModel(X, K=2)

        pd.DataFrame(Z, index=sample_names).to_csv(outfile+ "_ZIFA_Z.csv")
        pd.DataFrame(model_params['A'], index=np.concatenate(feature_names)).to_csv(outfile+ "_ZIFA_A.csv")
        pd.DataFrame(model_params['mus'], index=np.concatenate(feature_names)).to_csv(outfile+ "_ZIFA_mus.csv")
        pd.DataFrame(model_params['sigmas'], index=np.concatenate(feature_names)).to_csv(outfile+ "_ZIFA_sigmas.csv")
    else:
        from ZIFA import block_ZIFA
        Z, model_params = block_ZIFA.fitModel(X, K=2, p0_thresh=0.95)
        feature_names = np.array(feature_names)[(X == 0).sum(axis=0) / X.shape[0] <= 0.95]

        pd.DataFrame(Z, index=sample_names).to_csv(outfile+ "_ZIFA_Z.csv")
        pd.DataFrame(model_params['A'], index=np.concatenate(feature_names)).to_csv(outfile+ "_ZIFA_A.csv")
        pd.DataFrame(model_params['mus'], index=np.concatenate(feature_names)).to_csv(outfile+ "_ZIFA_mus.csv")
        pd.DataFrame(model_params['sigmas'], index=np.concatenate(feature_names)).to_csv(outfile+ "_ZIFA_sigmas.csv")


def train_PCA(input_data, feature_names, sample_names, outfile):
    from sklearn.decomposition import PCA
    
    X = []
    for m in range(len(data)):
        X.append(np.vstack(data[m]))  # concatenate samples across groups
    X = np.hstack(X) # concatenate features across views

    keep_sample = (~np.isnan(X)).sum(axis=1) > 0
    sample_names = np.concatenate(sample_names)[keep_sample]
    X = X[keep_sample, :]

    pca = PCA(n_components=2)
    pca.fit(X) 

    loadings = pd.DataFrame(pca.components_.transpose(), index=np.concatenate(feature_names))
    Z = pd.DataFrame(pca.transform(X), index=sample_names)

    loadings.to_csv(outfile + "_PCA_W.csv")
    Z.to_csv(outfile + "_PCA_Z.csv")


def train_CTF(input_data, feature_names, sample_names, times, outfile):
    
    from gemelli.ctf import ctf
    import biom
    # requires as input biom table and data frame with group and time annotations

    metadf = pd.DataFrame({'sample_id' : np.concatenate(sample_names)})
    featdf = pd.DataFrame({'feature_id' : np.concatenate(feature_names)})
    metadf[['group', 'month']] = metadf['sample_id'].str.split('_',2, expand = True)
    assert np.all(np.array(metadf.month, dtype = 'int') == np.concatenate(times).flatten())

    X = []
    for m in range(len(data)):
        X.append(np.vstack(data[m]))  # concatenate samples across groups
    X = np.hstack(X) # concatenate features across views

    # missing values in CTF encoded as zeros:
    X[np.isnan(X)] = 0

    table = biom.table.Table(X.transpose(), observation_ids = featdf.feature_id, sample_ids = metadf.sample_id)
    metadf.set_index('sample_id', inplace = True)
    featdf.set_index('feature_id', inplace = True)
    
    assert np.all(metadf.index.tolist() == table.ids("sample"))
    assert np.all(featdf.index.tolist() == table.ids("observation"))

    ord_res, state_ordn, dists, straj, ftraj = ctf(table, 
        sample_metadata = metadf,
        individual_id_column = 'group',
        state_column = 'month',
        n_components = 2,
        min_sample_count = 0,
        min_feature_count = 0,
        max_iterations_als = 5,
        max_iterations_rptm = 5,
        n_initializations = 5,
        feature_metadata = featdf)

    straj.to_csv(outfile + "_CTF_straj.csv")
    ftraj.to_csv(outfile + "_CTF_ftraj.csv")
    state_ordn.samples.to_csv(outfile + "_CTF_state_ord_samples.csv")
    state_ordn.features.to_csv(outfile + "_CTF_state_ord_features.csv")
    ord_res.samples.to_csv(outfile + "_CTF_subj_ord_samples.csv")
    ord_res.features.to_csv(outfile + "_CTF_subj_ord_features.csv")



if __name__ == '__main__':
    
    # I/O
    metadir = "data/processed_data/"
    
    # options
    method = sys.argv[1]
    Ndown = int(sys.argv[2])
    seed = int(sys.argv[3])
    preproc = sys.argv[4]


    # load data
    if preproc == "RCLR":
        assert method in ["MEFISTO", "MOFA", "CTF"] # requires handling of missing values and negative inputs
        if method == "CTF": # CTF does RCLR internally, start from corresponding count matrix
            datadir = "data/input_data_counts/"
            data, times, group_names, sample_names, feature_names, view_names = load_data_views(datadir, metadir)
        else:
            datadir = "data/input_data_rclr/"
            data, times, group_names, sample_names, feature_names, view_names = load_data_views(datadir, metadir)        
        outdir = "out/"
        center_groups = False
    else:
        print("Preprocessing not understood")
        sys.exit()

    # mask some samples at random
    if Ndown > 0:
        np.random.seed(seed)
        # only mask samples that have been observed in given data
        X = []
        for m in range(len(data)):
            X.append(np.vstack(data[m]))  # concatenate samples across groups
        X = np.hstack(X) # concatenate features across views
        keep_sample = (~np.isnan(X)).sum(axis=1) > 0
        samples2mask = np.concatenate(sample_names)[keep_sample]
        samples2mask = np.random.choice(samples2mask, Ndown, replace=False)
        for g in range(len(group_names)):
            for m in range(len(view_names)):
                data[m][g][np.isin(sample_names[g], samples2mask),:] = np.nan

    if method == "MOFA":
        train_MOFA(input_data = data,
            times = times,
            group_names = group_names,
            feature_names = feature_names,
            sample_names = sample_names,
            view_names = view_names,
            use_GP = False,
            center_groups = center_groups,
            outfile = outdir + "MOFA_%s_%s.hdf5" % (Ndown, seed))

    elif method == "MEFISTO":
        train_MOFA(input_data = data,
            times = times,
            group_names = group_names,
            feature_names = feature_names,
            sample_names = sample_names,
            view_names = view_names,
            use_GP = True,
            center_groups = center_groups,
            outfile = outdir + "MEFISTO_%s_%s.hdf5" % (Ndown, seed))
    
    elif method == "MEFISTO_nogroups":
        train_MOFA(input_data = data,
            times = times,
            group_names = group_names,
            feature_names = feature_names,
            sample_names = sample_names,
            view_names = view_names,
            use_GP = True,
            center_groups = center_groups,
            model_groups = False,
            outfile = outdir + "MEFISTO_nogroups_%s_%s.hdf5" % (Ndown, seed))

    elif method == "PCA":
        train_PCA(input_data = data,
            feature_names = feature_names,
            sample_names = sample_names,
            outfile = outdir + "PCA_%s_%s" % (Ndown, seed))

    elif method == "ZIFA":
        train_ZIFA(input_data = data,
            feature_names = feature_names,
            sample_names = sample_names,
            outfile = outdir + "ZIFA_%s_%s" % (Ndown, seed))
        
    elif method == "NMF":
         train_NMF(input_data = data,
                    feature_names = feature_names,
                    sample_names = sample_names,
                    outfile = outdir + "NMF_%s_%s" % (Ndown, seed))

    elif method == "CTF":
        train_CTF(input_data = data,
                feature_names = feature_names,
                sample_names = sample_names,
                times = times,
                outfile = outdir + "CTF_%s_%s" % (Ndown, seed))

    else:
        print("Method not implemented.")
        sys.exit()

