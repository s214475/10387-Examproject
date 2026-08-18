#!/bin/bash
#BSUB -J abtem_multislice
#BSUB -q hpc
#BSUB -W 06:00
#BSUB -n 16
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000MB]"
#BSUB -o /zhome/f1/0/167759/10387/data/logs/abtem_%J.out
#BSUB -e /zhome/f1/0/167759/10387/data/logs/abtem_%J.err

cd /zhome/f1/0/167759/10387
source .venv/bin/activate

python /zhome/f1/0/167759/10387/python/multislice.py