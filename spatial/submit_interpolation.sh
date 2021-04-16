. /software/anaconda3/2019.07/etc/profile.d/conda.sh
conda activate mefisto


for seed in 13421 14390 93243 434432 12334 58324 1244 18484 99732 41741
do
	bsub -n 1 -R "rusage[mem=1500MB]" -W 10:00 -e e_interpol_500 -o o_interpol_500 python interpolation_grid.py $seed 250 $frac_features 500
	bsub -n 1 -R "rusage[mem=2000MB]" -W 15:00 -e e_interpol_1000 -o o_interpol_1000 python interpolation_grid.py $seed 250 $frac_features 1000
	bsub -n 1 -R "rusage[mem=4000MB]" -W 20:00 -e e_interpol_full -o o_interpol_full python interpolation_grid.py $seed 250 $frac_features 5000
done


