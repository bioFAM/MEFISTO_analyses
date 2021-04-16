# Script preprocesses the ECAM data starting from the following files downloaded from https://codeocean.com/capsule/6494482/tree/v1:
# * samples metadata: /data/ECAM-Qiita-10249/10249_20180418-081211.txt
# * count table: /data/ECAM-Qiita-10249/table.biom
# * feature metadata: /data/ECAM-Qiita-10249/q2-analysis/taxonomy.qza
# These files are stored in the relative path data/raw_data.

# The preprocessing requires qiime (see setup_qiime.sh) and is done in a separate conda environemnt: conda activate qiime2-2020.8

import biom
import numpy as np
import pandas as pd
from biom import Table
from scipy.spatial import distance
import qiime2 as q2
import os
from gemelli.preprocessing import rclr

data_dir = 'data/raw_data/'

# load data
counts = biom.load_table(data_dir + 'table.biom')
samples_meta = pd.read_csv(data_dir + '10249_20180418-081211.txt',sep='\t',index_col=0)
features_meta = q2.Artifact.load(data_dir + 'taxonomy.qza')
features_meta = features_meta.view(q2.Metadata).to_dataframe()


# filter samples to children samples of type Stool_Stabilizer within the first 24 months of life 
samples_meta.index = samples_meta.index.astype(str)
samples_meta = samples_meta.replace('na',np.nan)
samples_meta = samples_meta.dropna(subset=['host_subject_id','month', 'sampletype', 'mom_child'])
samples_meta = samples_meta[samples_meta.sampletype == 'Stool_Stabilizer']
samples_meta['month'] = [int(x) for x in samples_meta.month]
samples_meta = samples_meta[samples_meta['month'].between(0, 24)] # note: CTF drops also 13, 17 and 21 (and posthoc (after the decomposition) 6,15,19)
samples_meta = samples_meta[samples_meta.mom_child == "C"] # note: CTF use C and M samples

# months per child
months_per_child = {k : sorted(set(df['month'])) for k, df in samples_meta.groupby('host_subject_id')}
n_months_per_child = [len(v) for k,v in months_per_child.items()]
print("Number of months per child: ", n_months_per_child)
print("Number of children: ", len(n_months_per_child))
children2keep = [k for k,v in months_per_child.items() if len(v) > 0] # keep all children with at least one measurement
samples_meta = samples_meta[samples_meta.host_subject_id.isin(children2keep)]



# filter counts to the retained samples 
counts.filter(samples_meta.index, axis = "sample", inplace = True)

# filter taxa without any counts
taxa2keep = counts.ids('observation')[counts.sum('observation') > 0]
# drop and filter
counts.filter(taxa2keep, axis='observation', inplace=True)

print("Keeping", counts.shape[1], "samples and", counts.shape[0], "features.")

# turn to data.frame (features x samples)
counts_df = pd.DataFrame(counts.matrix_data.toarray(),
                    counts.ids('observation'),
                    counts.ids('sample'))


# match indices for counts table and sample/feature meta 
samples_meta = samples_meta.reindex(counts_df.columns)
samples_meta.reset_index(inplace = True)
samples_meta.rename(columns={'index': 'SampleID'}, inplace = True)
samples_meta[['study_id','child_id','time_id']] = samples_meta['SampleID'].str.split(".",expand=True)

# subset feature table and match index
features_meta = features_meta.reindex(counts_df.index)

# save processed data
if not os.path.exists('data/processed_data'):
    os.makedirs('data/processed_data')
samples_meta.to_csv('data/processed_data/samples_metadata.csv', header = samples_meta.keys().get_values())
counts_df.index.name = 'SampleID'
counts_df.to_csv('data/processed_data/counts_matrix.csv')
features_meta.to_csv('data/processed_data/features_metadata.csv')

# print some data summaries
print("Number of samples per delivery mode:\n", samples_meta.delivery.value_counts())
print("Months with at least one sample:", set(samples_meta.month))
print("Number of samples by deliver mode and month:\n")
print(samples_meta.groupby(['month', 'delivery']).size())

# correct for library size and log transform
# ls = counts_df.sum(axis = 0)
# norm_df = np.log(counts_df / ls[None,:] * 10000 + 1)
# norm_df.to_csv('data/processed_data/log_normalized_matrix.csv')

# normalize data by RCLR (Martino et al)
rclr_mat = rclr(np.asarray(counts_df))
rclr_df = pd.DataFrame(rclr_mat, counts_df.index, counts_df.columns)
rclr_df.to_csv('data/processed_data/rclr_mat.csv')

# filter
# roughly 98% of taxas are zero per sample, keep only features present in at least 5 samples
n_min_samples = 5
counts_df_filtered = counts_df[(counts_df!=0).sum(axis = 1) >= n_min_samples]
rclr_df_filtered = rclr_df[(counts_df!=0).sum(axis = 1) >= n_min_samples]
# norm_df_filtered = norm_df[(counts_df!=0).sum(axis = 1) >= n_min_samples]

print("After filtering to minimum of", n_min_samples, "samples", counts_df_filtered.shape[0], "features are kept.")
print("Minimal number of non-zero species in a sample:", (counts_df_filtered > 0).sum(axis = 0).min())

counts_df_filtered.to_csv('data/processed_data/counts_matrix_filtered.csv')
# norm_df_filtered.to_csv('data/processed_data/log_normalized_matrix_filtered.csv')
rclr_df_filtered.to_csv('data/processed_data/rclr_mat_filtered.csv')



