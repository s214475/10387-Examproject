#!/bin/bash
#BSUB -J abtem_multislice
#BSUB -q hpc
#BSUB -W 08:00
#BSUB -n 32
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000MB]"
#BSUB -o /zhome/f1/0/167759/10387/data/logs/pt_wave_%J.out
#BSUB -e /zhome/f1/0/167759/10387/data/logs/pt_wave_%J.err

cd /zhome/f1/0/167759/10387
source .venv/bin/activate

python /zhome/f1/0/167759/10387/python/multislice.py