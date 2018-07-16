#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed 
#SBATCH --ntasks=1
#SBATCH --time=30:00:00
#SBATCH --mem=64GB
#SBATCH --output=/dors/capra_lab/users/fongsl/broadly_active_enhancers/data/

#SBATCH --job-name=20180716_analyze_hg38_100way_x_roadmap_each_tissue.py
#SBATCH --account=capra_lab
#SBATCH --error=%J.stdout
#SBATCH --output=%J.stderr

ml load Anaconda3
ml load  GCC/5.4.0-2.26  BEDTools/2.26.0

python /dors/capra_lab/users/fongsl/broadly_active_enhancers/bae_git/roadmap_each/20180716_analyze_hg38_100way_x_roadmap_each_tissue.py