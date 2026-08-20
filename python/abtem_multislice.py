from pathlib import Path
import os
import time
import logging

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("FFTW_NUM_THREADS", "1")

import abtem
from ase.io import read

base = Path("/zhome/f1/0/167759/10387/data")
model_dir = base / "models"
exitwave_dir = base / "exit_waves"
log_dir = base / "logs"

for folder in [model_dir, exitwave_dir, log_dir]:
    folder.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    filename=log_dir / "abtem_multislice.log",
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

n_workers = int(os.environ.get("LSB_DJOB_NUMPROC", "4"))

abtem.config.set({
    "dask.scheduler": "threads",
    "mkl.threads": 1,
    "fftw.threads": 1,
})

logging.info(f"Using {n_workers} workers")
logging.info(f"Running on host: {os.uname().nodename}")

lc = 3.92  # replace with your actual value if different

atom = read(model_dir / "Pt5nm_cov.traj")

energy = 300e3
slice_thickness = lc / 8
sampling = lc / 256

plane_wave = abtem.PlaneWave(energy=energy)

logging.info("Creating potential")
atom_pot = abtem.Potential(
    atom,
    slice_thickness=slice_thickness,
    sampling=sampling,
)

logging.info("Starting multislice")
t0 = time.perf_counter()

atom_ewave = plane_wave.multislice(atom_pot).compute(
    scheduler="threads",
    num_workers=n_workers,
)

elapsed = time.perf_counter() - t0
logging.info(f"Finished multislice in {elapsed / 60:.2f} min")

output = exitwave_dir / "Pt5nm_ewave.zarr"
atom_ewave.to_zarr(str(output), overwrite=True)

logging.info(f"Saved exit wave to {output}")
print(f"Finished in {elapsed / 60:.2f} min")
print(f"Saved to {output}")