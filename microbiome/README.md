# MEFISTO: Application to longitudinal microbiome data


In this folder you can find all script for the microbiome application of MEFISTO.

In particular,

* `preprocess_data.py` is used to preprocess the microbiome data for application of MEFISTO. This is based on the processing performed in [Martino et al](https://www.nature.com/articles/s41587-020-0660-7). As this requires `qiime2` you can use `setup_qiime.sh` to install the required dependencies.
* `run_ECAM.py` is the main script to prepare and train MEFISTO or alternative methods on the data and is called in `submit_training.sh`.
* `microbiome_analysis.Rmd` performs the downstream analysis of the trained model.
* `eval_stability.Rmd` evaluates the factor stability of MEFISTO to masking of some time points in comparison to other methods
* `utils.R` and `utils_stability.R` contain some helper function


The preprocessing is based on the code provided by [Martino et al](https://www.nature.com/articles/s41587-020-0660-7). Required are the following files downloaded from https://codeocean.com/capsule/6494482/tree/v1:

* samples metadata: /data/ECAM-Qiita-10249/10249_20180418-081211.txt
* count table: /data/ECAM-Qiita-10249/table.biom
* feature metadata: /data/ECAM-Qiita-10249/q2-analysis/taxonomy.qza

These files are stored in the relative path data/raw_data.