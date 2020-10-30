. /software/anaconda3/2019.07/etc/profile.d/conda.sh
conda activate mefisto
bsub -n 1 -R "rusage[mem=400MB]" -W 00:10 -e e_toy -o o_toy python toy_data.py time_aware
bsub -n 1 -R "rusage[mem=400MB]" -W 00:05 -e e_toy -o o_toy python toy_data.py time_nonaware
