. /software/anaconda3/2019.07/etc/profile.d/conda.sh
conda activate mefisto


#experiments=(noise G N missing align sparse) 
experiments=(align) 
for type in ${experiments[@]}
do
   for seed in 41741 13421 14390 93243 434432 12334 58334 1244 18484 99732
   do
       if [ $type = N ] 
       then
          for method in MOFA2 MEFISTO 
          do
            bsub -n 1 -R "rusage[mem=2500MB]" -W 60:00 -e e_$type -o o_$type python grid.py $type $seed $method
          done
       elif [ $type = G ] 
       then
          for method in MOFA2 MEFISTO
          do
            bsub -n 1 -R "rusage[mem=2500MB]" -W 24:00 -e e_$type -o o_$type python grid.py $type $seed $method
          done
       elif [ $type = align ]
       then
          for method in MOFA2 MEFISTO MEFISTO+align
          do
            bsub -n 1 -R "rusage[mem=2500MB]" -W 24:00 -e e_$type -o o_$type python grid.py $type $seed $method
          done
       elif [ $type = missing ] || [ $type = noise ]
       then
          for method in MOFA2 MEFISTO
          do
            bsub -n 1 -R "rusage[mem=1000MB]" -W 06:00 -e e_$type -o o_$type python grid.py $type $seed $method
          done
       elif [ $type = sparse ]
       then
          for method in MOFA2 MEFISTO MEFISTO_sparse
          do
            bsub -n 1 -R "rusage[mem=5000MB]" -W 60:00 -e e_$type -o o_$type python grid.py $type $seed $method
          done
        else
          echo "Argument type not understood"
       fi
   done
done
