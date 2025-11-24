#!/bin/bash
#SBATCH --job-name=antibody_test
#SBATCH --partition=gpu3090
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# load modules
module purge                                    
module load anaconda3
module load gromacs/2023.2-gcc-9.5.0-jzxesel            

# env with pandas numpy yaml biopython (peptide mode) and modeller
source activate ANUBI

# run
python ANUBI_main.py -i infile.yaml

