#!/bin/bash

#SBATCH --job-name=panoptes
#SBATCH --partition=gpu8_long,gpu4_long
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --begin=now
#SBATCH --time=8-00:00:00
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --mail-user=shini.chen@nyulangone.org
#SBATCH --output=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/logs/%j.out
#SBATCH --error=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/logs/%j.error
#SBATCH --gres=gpu:2

## to run:
# sbatch multi_train_test.sh X1


## activate env
source ~/.bashrc
conda activate panoptes_tf29

var=$1
echo "Running variant: $var"

#module purge
#module load cuda/11.8
#module load condaenvs/gpu/tensorflow2.9

cd /gpfs/data/proteomics/software/panoptes0/src/ # enter panoptes src folder to run train and test scripts

### train model
python train.py \
    --multi_gpu=True \
    --variant=$var \
    --split=True \
    --split_ratio=0.8,0.1,0.1 \
    --out_dir=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/results/$var \
    --tile_idx_dir=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/all_tiles_1qgain_filtered.csv \
    --label_df_dir=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/combined_HNSC_Arm1qGain_labels.csv \
    --lab_col='Arm_1q_Gain' \
    --max_epoch=20 \
    --batch_size=16 \
    --dropout=0.5 \
    --aux=True \
    --aux_weight=0.3 \
    --patience=8 \
    --seed=123456

### test model
#python test.py \
#    --multi_gpu=False \
#    --variant=$var \
#    --aux=True \
#    --out_dir=/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/results/$var/pred \
#    --tst_df=/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/results/$var/data/tst_tile_idx.csv \
#    --lab_level=Patient_ID \
#    --agg_level=Patient_ID \
#    --saved_model_dir=/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/results/$var/model/panoptes_weights_final.h5 \
#    --manifold_sample=20000 \
#    --seed=123456
