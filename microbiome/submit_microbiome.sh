. /software/anaconda3/2019.07/etc/profile.d/conda.sh
conda activate mefisto

bsub -n 1 -R "rusage[mem=1000MB]" -W 120:00 -e e_ECAM -o o_ECAM python run_ECAM.py
