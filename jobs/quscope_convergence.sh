#!/bin/bash

#BSUB -J quscope_conv
#BSUB -q hpc
#BSUB -n 4
#BSUB -W 720
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16000]"
#BSUB -o data/logs/quscope_conv_%J.out
#BSUB -e data/logs/quscope_conv_%J.err

# ============================================================
# Environment
# ============================================================

module load python3/3.11.13

source ~/10387/.venv/bin/activate

# ============================================================
# Threading
# ============================================================

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

# ============================================================
# Job information
# ============================================================

echo "=================================================="
echo "QuScope 3-slice convergence test"
echo "=================================================="

echo "Host:"
hostname

echo "Date:"
date

echo "LSF CPUs:"
echo "$LSB_DJOB_NUMPROC"

echo "Python:"
which python
python --version

echo "=================================================="

# ============================================================
# Run
# ============================================================

python ~/10387/python/quscope_convergence.py

status=$?

echo "=================================================="
echo "Finished with status: $status"
echo "Date:"
date
echo "=================================================="

exit $status
