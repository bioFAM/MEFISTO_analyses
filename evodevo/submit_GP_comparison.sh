. /software/anaconda3/2019.07/etc/profile.d/conda.sh
conda activate mefisto

for method in FA univGPs
do
  for seed in 93243 434432 12334 58324 1244 18484 99732 41741 13421 14390
  do
      for N in `seq 2 2 10`
      do
        for frac_features in 0.25 0.5 0.75 1
        do
             if [ $method = univGPs ]
             then
                bsub -n 1 -R "rusage[mem=700MB]" -W 06:00 -e e_interpol_comp_GP -o o_interpol_comp_GP python interpolation_grid_vs_GP.py $seed $N $method $frac_features
             else
                sleep 2s
                bsub -n 1 -R "rusage[mem=200MB]" -W 00:05 -e e_interpol_comp -o o_interpol_comp python interpolation_grid_vs_GP.py $seed $N $method $frac_features
             fi
        done
      done
  done
done