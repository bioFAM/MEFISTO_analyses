. /software/anaconda3/2019.07/etc/profile.d/conda.sh
conda activate mefisto

for K in 4 # number of factors
do
   for N in 500 1000 full # number of inducing points
   do
       if [ $N = 500 ]
       then
           bsub -n 1 -R "rusage[mem=1500MB]" -W 05:00 -e e_$N -o o_$N python run_brainST.py $N $K
       elif [ $N = 1000 ]
       then
           bsub -n 1 -R "rusage[mem=2000MB]" -W 10:00 -e e_$N -o o_$N python run_brainST.py $N $K
        else
           bsub -n 1 -R "rusage[mem=4000MB]" -W 15:00 -e e_$N -o o_$N python run_brainST.py $N $K
       fi
   done
done
