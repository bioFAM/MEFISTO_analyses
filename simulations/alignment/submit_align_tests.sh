#python alignment_tests.py 41741 imbalance

. /software/anaconda3/2019.07/etc/profile.d/conda.sh
conda activate mofa2

for seed in 13421 14390 93243 434432 12334 58324 1244 18484 99732 41741
do
    bsub -n 1 -R "rusage[mem=300MB]" -W 00:20 -e e_align_%J -o o_align_%J python alignment_tests.py $seed nonsmooth
    bsub -n 1 -R "rusage[mem=2000MB]" -W 00:20 -e e_align_%J -o o_align_%J python alignment_tests.py $seed imbalance
    bsub -n 1 -R "rusage[mem=2000MB]" -W 00:20 -e e_align_%J -o o_align_%J python alignment_tests.py $seed imbalance_weighted
#    bsub -n 1 -R "rusage[mem=2000MB]" -W 00:40 -e e_align_%J -o o_align_%J python alignment_tests.py $seed imbalance_extra
#    bsub -n 1 -R "rusage[mem=2000MB]" -W 00:40 -e e_align_%J -o o_align_%J python alignment_tests.py $seed imbalance_extra_weighted

done