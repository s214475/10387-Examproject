---
# 10387-Examproject
--- 

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

> .venv/

The environment is activated with:

> source ~/10387/.venv/bin/activate

To be allowed to use uv do

> export PATH="$PWD/.uv-tool/bin:$PATH"
