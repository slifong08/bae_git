#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed 
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=24GB
#SBATCH --output=/dors/capra_lab/users/fongsl/broadly_active_enhancers/data/

ml load Anaconda3
ml load GCC/5.4.0-2.26
ml load Intel/2016.3.210-GCC-5.4.0-2.26
ml load PHAST/1.4

python /dors/capra_lab/users/fongsl/bae_git/alignments/maf_alignment_pipeline_data.py