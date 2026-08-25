#!/bin/bash

#SBATCH --job-name=HNSC_Train
#SBATCH --partition=gpu4_medium
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --begin=now
#SBATCH --time=3-00:00:00
#SBATCH --mem=64GB
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END
#SBATCH --mail-user=shini.chen@nyulangone.org
#SBATCH --output=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/HNSC.HPVneg_9p_loss/logs/%j.out
#SBATCH --error=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/HNSC.HPVneg_9p_loss/logs/%j.error

## usage:
# sbatch multi_train.sh X1 CPTAC
# sbatch multi_train.sh X1 TCGA

# --- Activate environment ---
source ~/.bashrc
conda activate panoptes_tf29

# --- Load CUDA ---
module purge
module load cuda/11.8
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

# --- Parse input ---
variant=$1
dataset=$2
echo "Running variant: $variant on dataset: $dataset"

# --- Set dataset-specific paths ---
if [ "$dataset" == "CPTAC" ]; then
    label_file="/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/HNSC.HPVneg_9p_loss/cptac_HNSC_Arm9pLoss_labels.csv"
    tile_file="/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/HNSC.HPVneg_9p_loss/all_tiles_9ploss_filtered_CPTAC.csv"
elif [ "$dataset" == "TCGA" ]; then
    label_file="/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/HNSC.HPVneg_9p_loss/tcga_HNSC_Arm9pLoss_labels.csv"
    tile_file="/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/HNSC.HPVneg_9p_loss/all_tiles_9ploss_filtered_TCGA.csv"
else
    echo "Error: dataset must be either CPTAC or TCGA"
    exit 1
fi

# --- Change to panoptes code directory ---
cd /gpfs/data/proteomics/software/panoptes0/src/

# --- Train model ---
python train.py \
    --multi_gpu=False \
    --variant=$variant \
    --split=True \
    --split_ratio=0.8,0.1,0.1 \
    --out_dir=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/HNSC.HPVneg_9p_loss/results/${dataset}_$variant \
    --tile_idx_dir=$tile_file \
    --label_df_dir=$label_file \
    --lab_col='Arm_9p_Loss' \
    --max_epoch=20 \
    --batch_size=16 \
    --dropout=0.5 \
    --aux=True \
    --aux_weight=0.3 \
    --patience=8 \
    --seed=123456