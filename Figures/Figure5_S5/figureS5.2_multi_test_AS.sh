#!/bin/bash

#SBATCH --job-name=HNSC 
#SBATCH --partition=gpu8_short,gpu4_short
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --begin=now
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --mail-user=shini.chen@nyulangone.org
#SBATCH --output=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/HNSC/logs/%j.out
#SBATCH --error=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/HNSC/logs/%j.error
#SBATCH --gres=gpu:1

### Activate environment
source ~/.bashrc
conda activate panoptes_tf29

### Load CUDA
module purge
module load cuda/11.8
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

### Set working directory
cd /gpfs/data/proteomics/software/panoptes0/src/

### Run test
variant=X1
echo "Running variant: $variant"


# Lower perplexity in test.py:
python test.py \
    --multi_gpu=False \
    --variant=$variant \
    --aux=True \
    --out_dir='/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/HNSC/pred' \
    --tst_df=/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/results/X1/data/tst_tile_idx.csv  \
    --label_df_dir='/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/HNSC_HPVneg_globalAS_labels.csv' \
    --lab_col='High_Global_AS' \
    --lab_level='Patient_ID' \
    --agg_level='Patient_ID' \
    --saved_model_dir='/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/results/X1/model/panoptes_weights_best.h5' \
    --manifold_sample=0 \
    --seed=123456
