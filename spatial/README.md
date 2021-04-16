# MEFISTO: Application to spatial transcriptomics data


In this folder you can find all scripts for the application of MEFISTO to spatial transcriptomics data.

In particular,

* `run_brainST.py` is the main script to prepare and train MEFISTO on the data, `submit_spatial.sh` can be used to run all required configurations to reproduce the figures in the manuscript. 
* `spatial_transcriptomics.Rmd` prepares the input data for training with `run_brainST.py` and once trained performs the downstream analyses.
* `interpolation_grid.py`, `submit_interpolation.sh`, and `eval_interpol.Rmd` are used to evaluate the interpolation accuracy of the model for spatial locations which are missing.

The data is taken from the SeuratData package ('http://seurat.nygenome.org/src/contrib/stxBrain.SeuratData_0.1.1.tar.gz').

