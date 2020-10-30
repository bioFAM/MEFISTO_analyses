. /software/anaconda3/2019.07/etc/profile.d/conda.sh
conda activate mefisto
bsub -n 1 -R "rusage[mem=500MB]" -W 00:30 -e e_align_example -o o_align_example python grid.py align_example 13421 MEFISTO+align
