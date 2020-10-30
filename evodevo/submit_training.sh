. /software/anaconda3/2019.07/etc/profile.d/conda.sh
conda activate mefisto

# number of tissue-time combinations to mask at random
for Ndown in 0 1 3 5
do
	if [ $Ndown = 0 ] 
	then
		# only one seed  - all produce the same results as no subsampling
		bsub -n 1 -R "rusage[mem=300MB]" -W 01:00 -e e_downsample_all -o o_downsample_all python train_evodevo.py $Ndown 0 all 5
	else
		for seed in 345987 5328401 30481
   		do
      		bsub -n 1 -R "rusage[mem=300MB]" -W 01:00 -e e_downsample_all -o o_downsample_all python train_evodevo.py $Ndown $seed all 5
   		done
	fi
done
