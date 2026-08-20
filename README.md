# 10387-Examproject

This repository contains my mini project for DTU course **10387**.
Please have mercy on me, only look in jupyter notebook Notebooks/tutorial.ipynb

---
## How to use DTU HPC Cluster
---
The simulations are designed to run on DTU's HPC infrastructure using the LSF batch system.

The project uses a Python virtual environment located at:

```bash
cd ~/10387
```
The environment is activated with:

```bash
source .venv/bin/activate
```

To be allowed to use uv do

```bash
export PATH="$PWD/.uv-tool/bin:$PATH"
```

To run the .py files use the shell scripts .sh 

```bash
bsub < jobs/[name_of_file].sh
```

To see how the job is going use

```bash
bjobs
```