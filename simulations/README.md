# MEFISTO: Application to simulated data


In this folder you can find all script for the simulation studies on MEFISTO.

In particular,

* `grid.py` is the main script to simulate data and train MEFISTO on the data for various different parameters, e.g. varying the number of groups, timepoints, missing samples, noise level.
* `submit_grid.sh` is used to call the `grid.py` function on all parameter configurations included in the manuscript
* `evaluate_simulations.Rmd` visualizes the results from the simulation studies by visualizing factor revocery, imputation MSE and hyperparameter recovery in the different experiments
* `make_align_example.sh` calls the function defined in `grid.py` to generate an illustrative example of the alignment functionality of MEFISTO
* `example_warping.Rmd` visualizes the results from this illustrative example model
* `toy_data.py`, `submit_toy.sh` and `toy_data.Rmd` are used to simulate and evaluate an illustrative example to compare time-aware MEFISTO and time-agnostic MOFA on a single data set (as in Figure 1)
* `utils_simulation.R` contains some helper functions
* `alignment/` contains the files for additional tests on the alignment
* `univGP/` contains the files for additional comparisons of univariate GP regression models in terms of interpolation accuracy

