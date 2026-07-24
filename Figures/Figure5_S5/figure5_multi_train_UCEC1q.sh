#!/bin/bash

#SBATCH --job-name=UCEC_1q
#SBATCH --partition=gpu4_medium
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --begin=now
#SBATCH --time=3-00:00:00
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --mail-user=shini.chen@nyulangone.org
#SBATCH --output=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/UCEC_1q_gain/logs/%j.out
#SBATCH --error=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/UCEC_1q_gain/logs/%j.error
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4

## activate env
source ~/.bashrc
conda activate panoptes_tf29

# Load CUDA 11.8 module
module purge
module load cuda/11.8

# Optional: fix LD_LIBRARY_PATH if still missing CUDA libs
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

var=$1
echo "Running variant: $var"

cd /gpfs/data/proteomics/software/panoptes0/src/ # enter panoptes src folder to run train and test scripts

### train model
python train.py \
    --multi_gpu=False \
    --variant=$var \
    --split=True \
    --split_ratio=0.8,0.1,0.1 \
    --out_dir=/gpfs/data/proteomics/projects/Lisa/Image_Models/UCEC_MSI/results/$var \
    --tile_idx_dir=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/UCEC_1q_gain/tiles.csv \
    --label_df_dir=//gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/UCEC_1q_gain/arm1q_labels.csv \
    --lab_col='Arm_1q_Gain' \
    --max_epoch=20 \
    --batch_size=16 \
    --dropout=0.5 \
    --aux=True \
    --aux_weight=0.3 \
    --patience=8 \
    --seed=123456
