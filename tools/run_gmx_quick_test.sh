#!/bin/bash
# Quick GROMACS GPU timing test script
# Usage: ./run_gmx_quick_test.sh /path/to/file.tpr [ntmpi] [ntomp] [nsteps] [gpu_flags]

TPR=${1:-system.tpr}
NTMPI=${2:-1}
NTOMP=${3:-4}
NSTEPS=${4:-100}
GMX=/home/albert/apps/miniforge3/envs/gromacs/bin.AVX2_256/gmx
GPU_FLAGS=${5:-"-nb gpu -pme gpu -bonded gpu -update gpu"}

if [ ! -f "$TPR" ]; then
  echo "TPR file not found: $TPR"
  exit 1
fi

export OMP_NUM_THREADS=${NTOMP}
export MKL_NUM_THREADS=1

OUT=quick_gmx_test_$(basename "$TPR" .tpr).log
echo "Running quick GROMACS test on $TPR -> $OUT"

time $GMX mdrun -s "$TPR" -ntmpi $NTMPI -ntomp $NTOMP $GPU_FLAGS -nsteps $NSTEPS -v &> "$OUT"

echo "Done. Log saved to $OUT"
#!/bin/bash
# Quick GROMACS GPU timing test script
# Usage: ./run_gmx_quick_test.sh /path/to/file.tpr [ntmpi] [ntomp] [nsteps]

TPR=${1:-system.tpr}
NTMPI=${2:-1}
NTOMP=${3:-4}
NSTEPS=${4:-100}
GMX=/home/albert/apps/miniforge3/envs/gromacs/bin.AVX2_256/gmx
GPU_FLAGS=${5:-"-nb gpu -pme gpu -bonded gpu -update gpu"}

if [ ! -f "$TPR" ]; then
  echo "TPR file not found: $TPR"
  exit 1
fi

export OMP_NUM_THREADS=${NTOMP}
export MKL_NUM_THREADS=1

OUT=quick_gmx_test_$(basename "$TPR" .tpr).log
echo "Running quick GROMACS test on $TPR -> $OUT"

time $GMX mdrun -s "$TPR" -ntmpi $NTMPI -ntomp $NTOMP $GPU_FLAGS -nsteps $NSTEPS -v &> "$OUT"

echo "Done. Log saved to $OUT"
