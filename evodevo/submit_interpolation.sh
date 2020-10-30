. /software/anaconda3/2019.07/etc/profile.d/conda.sh
conda activate mefisto


for seed in 41741 13421 14390 93243 434432 12334 58324 1244 18484 99732
do
    for N in `seq 2 2 20`
    do
    	for Nviews in 3 4 5
    	do
    		bsub -n 1 -R "rusage[mem=300MB]" -W 01:00 -e e_interpol -o o_interpol python interpolation_grid.py $seed $N $Nviews
    	done
    done
done

