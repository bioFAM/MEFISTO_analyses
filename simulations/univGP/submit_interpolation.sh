. /software/anaconda3/2019.07/etc/profile.d/conda.sh
conda activate testGP

for seed in 12334 58334 1244 18484 99732 41741 13421 14390 93243 434432
do
  for method in MEFISTO univGPs MOFA2
  do
    for experiment in D smoothness missing N n_factors
    do
       if [ $method = univGPs ]
       then
           bsub -n 1 -R "rusage[mem=1000MB]" -W 05:00 -e e_interpol_$method -o o_interpol_$method python interpolation_grid.py $experiment $seed $method
       else
           bsub -n 1 -R "rusage[mem=500MB]" -W 00:10 -e e_interpol_$method -o o_interpol_$method python interpolation_grid.py $experiment $seed $method
       fi
    done
  done
done

