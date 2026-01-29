#!/bin/bash
#SBATCH --job-name=antibody_test
#SBATCH --partition=gpu3090
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# load modules
module purge                                    
module load anaconda3
# Prefer a modern GROMACS build with CUDA support (adjust module name/version as needed)
module load gromacs/2024.3-cuda

# Activate conda env used for this pipeline
source activate ANUBI

# Tune threads for GROMACS: one MPI rank per GPU and OpenMP threads per rank
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}
export MKL_NUM_THREADS=1

# Optional: pin threads for better performance (uncomment and adjust if supported)
# export GMX_PIN_THREADS=1

# Example gmx mdrun invocation (used internally by ANUBI). For manual runs:
# gmx mdrun -s system_Compl_MD.tpr -ntmpi 1 -ntomp ${OMP_NUM_THREADS} -nb gpu -pme gpu -bonded gpu -update gpu -pin on -pinoffset 0

# Run ANUBI pipeline
python ANUBI_main.py -i infile.yaml

