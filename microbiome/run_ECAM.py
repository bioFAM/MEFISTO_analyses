import pickle
import os
import numpy as np
import pandas as pd
import re
from mofapy2.run.entry_point import entry_point


seed = 2020
datadir= "data/processed_data/"

# Prepare input data for MOFA
# RCLR normalized data (mean over samples with multiple conditional overlaps)
# first dimension = samples; second dimension = features; [3..N] dimensions = conditions
rclr_counts = pickle.load(open(datadir + "ECAM_rclr_counts.p", "rb"))

# non-observed values have been set to 0
rclr_counts[rclr_counts == 0] = np.nan 
rclr_counts.shape

# reshape to list of groups
input_data = []
for i in range(rclr_counts.shape[0]):
    input_data.append(rclr_counts[i,:,:].transpose())
len(input_data)
input_data = [input_data] # list of 1 view


# load meta data
sample_meta = pd.read_csv(datadir + "sample_metadata.csv" )
feature_meta = pd.read_csv(datadir + "taxonomy.csv" )
table = pickle.load(open(datadir + "ECMA_table.p", "rb"))
tensor = pickle.load(open(datadir + "ECMA_tensor.p", "rb"))

# add feautre and view names:
tensor = pickle.load(open(datadir + "ECMA_tensor.p", "rb"))
feature_names = tensor.feature_order
view_names = ["microbiome"]

# add sample and group names:
sample_meta['baby_id'] = [re.split('\.',h)[1] for h in sample_meta['#SampleID']]
groupnames = [sample_meta[sample_meta.host_subject_id == subj].baby_id.values[0] for subj in tensor.subject_order]
sample_names = [[g + "-" + str(s) for s in tensor.condition_orders[0]] for g in groupnames ]

# subset as described in methods part and mother sample (?)
sample2remove = np.where([g == 'M050' for g in groupnames])[0].item()

# too few subjects at these months 
times2keep = np.where([t not in [6, 15, 19] for t in tensor.condition_orders[0]])[0]
del groupnames[sample2remove]
del sample_names[sample2remove]
del input_data[0][sample2remove]
sample_names = [np.array(s)[times2keep] for s in sample_names]
input_data = [[input_data[0][g][times2keep,:] for g in range(len(groupnames))]]

# covariate
months = np.asarray(tensor.condition_orders[0])[times2keep]
months = [months[:,None]] * len(groupnames)

# prepare MEFISTO model
ent = entry_point()
ent.set_data_options(center_groups=False)
ent.set_data_matrix(input_data, groups_names=groupnames,
                    features_names=[feature_names], samples_names=sample_names,
                    views_names=view_names)
ent.set_model_options(factors=2)
ent.set_train_options(seed=seed, convergence_mode = "medium")
ent.set_covariates(months, covariates_names="month")
ent.set_smooth_options(model_groups = True, warping=False, warping_ref=0, n_grid=10, opt_freq=50, start_opt = 50)

# Build and run the model
ent.build()
ent.run()

# interpolate
ent.predict_factor(new_covariates=ent.model.nodes["Sigma"].covariates)

# save
outfile = "out/ECAM_microbiome-seed_%s.hdf5" % seed
	
ent.save(outfile)