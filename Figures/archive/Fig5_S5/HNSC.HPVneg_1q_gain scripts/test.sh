#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --partition=gpu8_short,gpu4_short
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --begin=now
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --mail-user=shini.chen@nyulangone.org
#SBATCH --output=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/logs/%j.out
#SBATCH --error=/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/logs/%j.error
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
var=$1
echo "Running variant: $var"

python test.py \
    --multi_gpu=False \
    --variant=$var \
    --aux=True \
    --out_dir='/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/results/'$var'/pred' \
    --tst_df='/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/results/'$var'/data/tst_tile_idx.csv' \
    --lab_level='Patient_ID' \
    --agg_level='Patient_ID' \
    --saved_model_dir='/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/image_analysis/results/'$var'/model/panoptes_weights_final.h5' \
    --manifold_sample=20000 \
    --seed=123456