#!/bin/bash

#BSUB -J quscope_ms
#BSUB -q hpc
#BSUB -n 4
#BSUB -W 180
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16000]"
#BSUB -o data/logs/quscope_%J.out
#BSUB -e data/logs/quscope_%J.err


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
# Run
# ============================================================

echo "=================================================="
echo "QuScope multislice"
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

python ~/10387/python/quscope_multislice.py

status=$?

echo "=================================================="
echo "Finished with status: $status"
echo "Date:"
date
echo "=================================================="

exit $status