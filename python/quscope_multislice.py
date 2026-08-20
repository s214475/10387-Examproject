from pathlib import Path
import os
import time
import logging

import numpy as np

from quscope.quantum_ctem import (
    QuantumMultisliceCircuit,
    QuantumMultisliceParameters,
)


# ============================================================
# Paths
# ============================================================

BASE = Path("/zhome/f1/0/167759/10387/data")

POTENTIAL_FILE = (
    BASE / "models" / "Pt5nm_cov_quscope.npz"
)

OUTPUT_FILE = (
    BASE / "exit_waves" / "Pt5nm_quscope_result.npz"
)

LOG_DIR = BASE / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)

LOG_FILE = LOG_DIR / "quscope_multislice.log"


# ============================================================
# Thread configuration
# ============================================================

# Match the 4-core LSF job.
n_threads = int(os.environ.get("LSB_DJOB_NUMPROC", "4"))

os.environ.setdefault("OMP_NUM_THREADS", str(n_threads))
os.environ.setdefault("MKL_NUM_THREADS", str(n_threads))
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(n_threads))
os.environ.setdefault("NUMEXPR_NUM_THREADS", str(n_threads))


# ============================================================
# Logging
# ============================================================

logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

logging.info("Starting QuScope multislice calculation")
logging.info("Threads requested: %d", n_threads)
logging.info("Potential file: %s", POTENTIAL_FILE)
logging.info("Output file: %s", OUTPUT_FILE)


# ============================================================
# Load potential
# ============================================================

print("=" * 70)
print("QuScope quantum multislice")
print("=" * 70)

print(f"Loading potential:")
print(f"  {POTENTIAL_FILE}")

t0 = time.perf_counter()

data = np.load(POTENTIAL_FILE)

potentials = data["potentials"]

grid_size = int(data["grid_size"])
pixel_size = float(data["pixel_size"])
slice_thickness = float(data["slice_thickness"])

if "field_of_view" in data:
    field_of_view = float(data["field_of_view"])
else:
    field_of_view = grid_size * pixel_size

load_time = time.perf_counter() - t0


# ============================================================
# Validate potential
# ============================================================

n_slices = potentials.shape[0]

print()
print("Potential:")
print(f"  shape          = {potentials.shape}")
print(f"  number slices  = {n_slices}")
print(f"  grid size      = {grid_size}")
print(f"  pixel size     = {pixel_size:.6f} Å")
print(f"  slice thickness= {slice_thickness:.6f} Å")
print(f"  field of view  = {field_of_view:.3f} Å")

logging.info("Potential shape: %s", potentials.shape)
logging.info("Grid size: %d", grid_size)
logging.info("Pixel size: %.8f", pixel_size)
logging.info("Slice thickness: %.8f", slice_thickness)

if potentials.ndim != 3:
    raise ValueError(
        f"Expected potential with shape "
        f"(n_slices, N, N), got {potentials.shape}"
    )

if potentials.shape[1] != grid_size:
    raise ValueError(
        f"Potential x dimension ({potentials.shape[1]}) "
        f"does not match grid_size ({grid_size})"
    )

if potentials.shape[2] != grid_size:
    raise ValueError(
        f"Potential y dimension ({potentials.shape[2]}) "
        f"does not match grid_size ({grid_size})"
    )


# ============================================================
# Simulation parameters
# ============================================================

acceleration_voltage = 300e3
defocus = -500
cs = 1.3

print()
print("Simulation parameters:")
print(f"  acceleration voltage = {acceleration_voltage / 1e3:.0f} keV")
print(f"  grid size             = {grid_size}")
print(f"  pixel size            = {pixel_size:.6f} Å")
print(f"  slice thickness       = {slice_thickness:.6f} Å")
print(f"  defocus               = {defocus} Å")
print(f"  Cs                    = {cs}")


# ============================================================
# Build QuScope circuit
# ============================================================

print()
print("Building QuScope circuit...")

t0 = time.perf_counter()

params = QuantumMultisliceParameters(
    acceleration_voltage=acceleration_voltage,
    grid_size=grid_size,
    pixel_size=pixel_size,
    defocus=defocus,
    cs=cs,
    slice_thickness=slice_thickness,
)

circuit = QuantumMultisliceCircuit(params)

build_time = time.perf_counter() - t0

print(f"Circuit construction: {build_time:.2f} s")

logging.info(
    "Circuit construction time: %.3f s",
    build_time,
)


# ============================================================
# Run full multislice calculation
# ============================================================

print()
print("=" * 70)
print("STARTING FULL QUANTUM MULTISLICE")
print("=" * 70)

t0 = time.perf_counter()

result = circuit.simulate(potentials)

simulation_time = time.perf_counter() - t0

print()
print("=" * 70)
print("SIMULATION COMPLETE")
print("=" * 70)

print(
    f"Simulation time: "
    f"{simulation_time / 60:.2f} min"
)

logging.info(
    "Simulation time: %.3f s (%.2f min)",
    simulation_time,
    simulation_time / 60,
)


# ============================================================
# Inspect result
# ============================================================

print()
print("Result keys:")
print(result.keys())

for key, value in result.items():

    if isinstance(value, np.ndarray):

        print(
            f"  {key:15s}: "
            f"shape={value.shape}, "
            f"dtype={value.dtype}"
        )

    else:

        print(
            f"  {key:15s}: "
            f"{type(value)}"
        )


# ============================================================
# Save result
# ============================================================

OUTPUT_FILE.parent.mkdir(
    parents=True,
    exist_ok=True,
)

print()
print("Saving result...")

t0 = time.perf_counter()

np.savez_compressed(
    OUTPUT_FILE,

    # QuScope outputs
    statevector=result["statevector"],
    wave_function=result["wave_function"],
    amplitude=result["amplitude"],
    phase=result["phase"],
    intensity=result["intensity"],

    # Simulation metadata
    grid_size=grid_size,
    pixel_size=pixel_size,
    slice_thickness=slice_thickness,
    field_of_view=field_of_view,
    acceleration_voltage=acceleration_voltage,
    defocus=defocus,
    cs=cs,
    n_slices=n_slices,
)

save_time = time.perf_counter() - t0

print(f"Saved to:")
print(f"  {OUTPUT_FILE}")

print(f"Save time: {save_time:.2f} s")

logging.info(
    "Saved result to %s",
    OUTPUT_FILE,
)

logging.info(
    "Total simulation time: %.3f s",
    simulation_time,
)

print()
print("=" * 70)
print("DONE")
print("=" * 70)