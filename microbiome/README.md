# MEFISTO: Application to longitudinal microbiome data


In this folder you can find all script for the microbiome application of MEFISTO.

In particular,

* `preprocess_data.py` is used to preprocess the microbiome data for application of MEFISTO. This is based on the processing performed in [Martino et al](https://www.nature.com/articles/s41587-020-0660-7). AS this requires `qqime2` you can use `setup_qiime.sh` to install the required dependencies.
* `run_ECAM.py` is the main script to prepare and train MEFISTO on the data and is called in `submit_microbiome.sh`.
* `ECAM_analysis.Rmd` performs the downstream analysis of the trained model.

The preprocessing is based on the code provided by [Martino et al](https://www.nature.com/articles/s41587-020-0660-7) as a ‘Code Ocean’ capsule: https://doi.org/10.24433/CO.5938114.v1. , where you can also find the input data files (required are `table.qza`, `metadata-matched.tsv` and `taxonomy.qza`)