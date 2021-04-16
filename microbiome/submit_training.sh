. /software/anaconda3/2019.07/etc/profile.d/conda.sh
conda activate mefisto

 seedlist=(3428095 432521 4235 987642 4791 37180 123 49230 832950 13523)
 Nlist=(0 3 5 10 20 50 100 150)
for Ndown in ${Nlist[@]}
do
	if [ $Ndown = 0 ] 
	then
		# only one seed  - all produce the same results as no subsampling		
		bsub -n 1 -R "rusage[mem=500MB]" -W 00:15 -e e_MOFA -o o_MOFA python run_ECAM.py MOFA 0 0 RCLR
		bsub -n 1 -R "rusage[mem=1000MB]" -W 50:00 -e e_MEFISTO -o o_MEFISTO python run_ECAM.py MEFISTO 0 0 RCLR
		bsub -n 1 -R "rusage[mem=1000MB]" -W 01:00 -e e_CTF -o o_CTF python run_ECAM.py CTF 0 0 RCLR
	else
		for seed in ${seedlist[@]}
   		do
			bsub -n 1 -R "rusage[mem=500MB]" -W 00:15 -e e_MOFA -o o_MOFA python run_ECAM.py MOFA $Ndown $seed RCLR
			bsub -n 1 -R "rusage[mem=1000MB]" -W 50:00 -e e_MEFISTO -o o_MEFISTO python run_ECAM.py MEFISTO $Ndown $seed RCLR
			bsub -n 1 -R "rusage[mem=1000MB]" -W 01:00 -e e_CTF -o o_CTF python run_ECAM.py CTF $Ndown $seed RCLR
   		done
	fi
done
