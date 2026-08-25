#!/bin/bash

#SBATCH --job-name=panoptes_test
#SBATCH --partition=gpu8_short,gpu4_short
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --begin=now
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --mail-user=lizabeth.katsnelson@nyulangone.org
#SBATCH --output=test_%j.out
#SBATCH --error=test_%j.error
#SBATCH --gres=gpu:1

## to run:
# sbatch multi_train_test.sh X1

var=$1

module load condaenvs/gpu/tensorflow2.5
module load cuda/11.8

cd /gpfs/data/proteomics/software/panoptes0/src/ # enter panoptes src folder to run train and test scripts

### test model
python test.py \
    --multi_gpu=False \
    --variant=$var \
    --aux=True \
    --out_dir=/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/results/$var/pred \
    --tst_df=/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/results/$var/data/tst_tile_idx.csv \
    --lab_level=Patient_ID \
    --agg_level=Patient_ID \
    --saved_model_dir=/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/results/$var/model/panoptes_weights_final.h5 \
    --manifold_sample=20000 \
    --seed=123456
