# MEFISTO: Application to evodevo data


In this folder you can find all script for the evodevo application of MEFISTO.

In particular,

* `train_evodevo.py` is the main script to prepare and train MEFISTO on the full as well as on downsampled data, `submit_training.sh` can be used to run all required configurations to reproduce the figures in the manuscript.
* `evodevo_all.Rmd` prepares the input data for training analyses with `train_evodevo.py` and once trained performs the downstream analyses for application on mouse, rat, rabbit, human and opossum
* `interpolation_grid.py` is a script to mask view-timepoint combinations at random and evaluate the imputation/interpolation capacities of MOFA and MEFISTO. `submit_interpolation.sh` can be used  to run all required configurations to reproduce the figures in the manuscript.
* `eval_interpol.Rmd` visualizes the results from the interpolation experiments
* `utils_evodevo.R` contains some helper functions that are used in the Rmarkdowns

The data and supplementary information on the evodevo atlas is available from [Cardoso-Moreira et al, Nature 2019](https://www.nature.com/articles/s41586-019-1338-5). `preprocess_data.R` normalizes the count data from the original study for usage with MEFISTO.