#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed 
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=24GB
#SBATCH --output=/dors/capra_lab/users/fongsl/broadly_active_enhancers/data/

#SBATCH --job-name=get_hspecific_alignments
#SBATCH --account=capra_lab
#SBATCH --error=%J.stdout
#SBATCH --output=%J.stderr

ml load Anaconda3

python /dors/capra_lab/users/fongsl/broadly_active_enhancers/bae_git/alignments/get_hspecific_alignments_hg19.py 