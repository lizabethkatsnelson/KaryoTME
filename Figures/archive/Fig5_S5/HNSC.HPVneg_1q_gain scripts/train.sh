#!/bin/bash
#
#SBATCH --job-name=panoptes
#SBATCH --partition=gpu8_long,gpu4_long
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --begin=now
#SBATCH --time=8-00:00:00
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --mail-user=lizabeth.katsnelson@nyulangone.org
#SBATCH --output=hnsc_9p_loss_X1_%j.out
#SBATCH --error=hnsc_9p_loss_X1_%j.error
#SBATCH --gres=gpu:2

#module load condaenvs/gpu/tensorflow2.2
module load condaenvs/gpu/tensorflow2.5
module load cuda/11.8


cd /gpfs/data/proteomics/software/panoptes0/src/

python train.py \
    --multi_gpu=True \
    --variant=X1 \
    --split=True \
    --split_ratio=0.8,0.1,0.1 \
    --out_dir='/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/results/X1' \
    --tile_idx_dir='/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/all_tiles.csv' \
    --label_df_dir='/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/arm9p_loss_labels.csv' \
    --lab_col='Arm_9p_Loss' \
    --max_epoch=20 \
    --batch_size=16 \
    --dropout=0.7 \
    --base_model='InceptionResNetV2' \
    --feature_pool=False \
    --aux=True \
    --aux_weight=0.3 \
    --patience=4 \
    --seed=42
