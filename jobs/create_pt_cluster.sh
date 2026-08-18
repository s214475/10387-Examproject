#!/bin/bash
#BSUB -J create_pt_cluster
#BSUB -q hpc
#BSUB -W 01:00
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4000MB]"
#BSUB -o /zhome/f1/0/167759/10387/data/logs/create_pt_cluster_%J.out
#BSUB -e /zhome/f1/0/167759/10387/data/logs/create_pt_cluster_%J.err

cd /zhome/f1/0/167759/10387
source .venv/bin/activate

python /zhome/f1/0/167759/10387/python/create_pt_cluster.py