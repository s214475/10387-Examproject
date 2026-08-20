# 10387-Examproject

This repository contains my mini project for DTU course **10387**.
Please have mercy on me, only look in jupyter notebook Notebooks/tutorial.ipynb

---
## WorkFlow
---

10387/
│
├── data/
│   ├── models/
│   │   ├── Pt5nm_prime.traj
│   │   ├── Pt5nm_cov.traj
│   │   └── Pt5nm_cov_quscope.npz
│   │
│   ├── exit_waves/
│   │   └── ...
│   │
│   └── logs/
│       └── ...
│
├── jobs/
│   ├── create_pt_cluster.sh
│   ├── multislice.sh
│   ├── quscope_multislice.sh
│   └── test.sh
│
├── python/
│   ├── multislice.py
│   ├── quscope_model.py
│   ├── quscope_multislice.py
│   └── test.py
│
├── src/
│   └── pt_cluster/
│       └── functions.py
│
├── Notebooks/
│   ├── jungle.ipynb
│   └── tutorial.ipynb
│
├── preamp.txt
├── pyproject.toml
├── uv.lock
├── .gitignore
└── README.md
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