# 10387-Examproject

This repository contains my exam project for DTU course **10387**.
Please have mercy on me, only look in jupyter notebook the notebooks/tutorial.ipynb it is a jungle of thoughts.
The folder is used to verify codes/functions.
Too investigae the py files:  should be here inside 

---
## Download the repository
---
This project uses [git](https://git-scm.com/) to manedge files. First copy the repository to your working folder: 

```bash
git clone https://github.com/s214475/10387-Examproject
```

Then enter the repository 

```bash
cd 10387-Examproject
```

This project uses [uv](https://docs.astral.sh/uv/) to manedge the Python enviroment. When the repository is cloned the enviroment can be set using:

```bash
uv sync
```

---
## Repository structure 
---
```bash
Notebooks/
    jungle.ipynb                Do not open this, it is very weird oacked for data
    tutorial.ipynb              nly open this notebook, contains documentaion of the project. 

data/
    exit_waves/                 List of exit waves used for the analysis
    logs/                       Log of jobs with their error and output files of the job.
    models/                     Ase and quscope model of the Pt NP 

jobs/
    create_pt_cluster.sh        Shell file to create_pt_cluster.py
    abtem_multislice.sh         Shell file to initiate abtem_multislice.py
    quscope_multislice.sh       Shell file to initiate quscope_multislice.py

python/
    create_pt_cluster.py        Create an ase Pt NP with full CO coverage using bachelor thesis code.
    abtem_multislice.py         Calculate the clasic exit wave using abTEM and the ase model.
    quscope_convergence         Investigate runtime for different field of view resolutions.
    quscope_multislice.py       Calculate the exit wave using quscope model.
    quscope_model.py            Translate ase structure into quscope potential model.

src/                            Core Python source files and application logic.
    __pycache__                 Build artifacts and Python bytecode is automatically .gitignore
    __init__.py                 Used to define functions inside the folder as a package.
    functions.py                Add Bachelor thesis functions as a package for the env. 
```

---
## How to use DTU HPC Cluster
---
The simulations are designed to run on DTU's HPC infrastructure using the LSF batch system.

The project uses a Python virtual environment located at:

```bash
cd 10387-Examproject
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