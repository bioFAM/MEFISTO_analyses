. /software/anaconda3/2019.07/etc/profile.d/conda.sh
conda activate qiime2-2020.8

bsub -n 1 -R "rusage[mem=500MB]" -W 00:10 -e e_preprocess_ECAM_%J -o o_preprocess_ECAM_%J python preprocess_data.py

