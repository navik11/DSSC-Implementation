#!/bin/bash -l
#SBATCH --gres=gpu:1
#SBATCH --mem=20G
#SBATCH -J dssctest

module purge
conda activate /home/sachida/miniconda3/envs/dgl_env

# f=../data/sample_151507_anno.h5
# for i in {1..5}
#  do
# python -u run_DSSC.py --n_clusters -1 --data_file $f --save_dir out_151507 \
# --final_labels pred.csv --final_latent_file latent.csv --run $i \
# --ml_file ../example_constraints/sample_151507_mlFromMarks_test.txt --cl_file ../example_constraints/sample_151507_clFromMarks_test.txt
#  done

f=../data/Dataset2_osmFISH.h5ad
for i in {1..7}
 do
python -u run_DSSC2.py --n_clusters 11 --data_file $f --save_dir out_processed_osmFISH \
--select_genes 30 --encodeLayer 16 --decodeLayer 16 --z_dim 8 \
--final_labels pred.csv --final_latent_file latent.csv --run $i \
--ml_file ../example_constraints/sample_osmFish_mlFromMarks_test.txt --cl_file ../example_constraints/sample_osmFish_clFromMarks_test.txt
 done
