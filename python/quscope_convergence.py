"""
QuScope 3-slice transverse-sampling convergence test.

Runs the SAME three middle physical slices of the ASE Pt/C/O model at:

    256 x 256, 0.500 Å/pixel
    512 x 512, 0.250 Å/pixel
    1024 x 1024, 0.125 Å/pixel

The transverse field of view is kept fixed at 128 Å in all cases.
This preserves substantial vacuum around the ~56 Å particle.

The z slice thickness is kept equal to the abTEM calculation:

    slice_thickness = lc / 8

For each resolution:
    1. Build the 3-slice Kirkland potential from the ASE model.
    2. Run QuantumMultisliceCircuit.
    3. Save the result immediately to an NPZ file.

Results are saved separately so that completed resolutions are not lost
if a later, more expensive resolution is terminated.
"""

from pathlib import Path
from collections import Counter
import time
import os

import numpy as np
from ase.io import read

from quscope.ctem import KirklandPotential
from quscope.quantum_ctem import (
    QuantumMultisliceCircuit,
    QuantumMultisliceParameters,
)


# ============================================================
# Paths
# ============================================================

BASE = Path("/zhome/f1/0/167759/10387/data")

MODEL_FILE = BASE / "models" / "Pt5nm_cov.traj"
OUTPUT_DIR = BASE / "exit_waves" / "quscope_convergence"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Physical parameters
# ============================================================

LC = 3.92
ACCELERATION_VOLTAGE = 300e3
SLICE_THICKNESS = LC / 8       # ~0.49 Å
FIELD_OF_VIEW = 128.0          # Å

DEFOCUS = -500.0               # Å
CS = 1.3                       # mm

# Fixed FOV, increasing transverse sampling:
RESOLUTIONS = {
    256: FIELD_OF_VIEW / 256,   # 0.5 Å/pixel
    512: FIELD_OF_VIEW / 512,   # 0.25 Å/pixel
    1024: FIELD_OF_VIEW / 1024, # 0.125 Å/pixel
}

# Number of physical z-slices to simulate
N_TEST_SLICES = 3


# ============================================================
# Threading
# ============================================================

n_workers = int(os.environ.get("LSB_DJOB_NUMPROC", "4"))

os.environ.setdefault("OMP_NUM_THREADS", str(n_workers))
os.environ.setdefault("MKL_NUM_THREADS", str(n_workers))
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(n_workers))
os.environ.setdefault("NUMEXPR_NUM_THREADS", str(n_workers))


# ============================================================
# Kirkland parameters
# ============================================================

kirkland = KirklandPotential()

# Pt Kirkland parameterization used in the QuScope calculation.
# Source: E. J. Kirkland, Advanced Computing in Electron Microscopy,
# 2nd ed., Springer (2010), as provided through the abTEM parameter data.
kirkland.params_dict["Pt"] = [
    [0.98469794, 2.73987079, 3.61696715],       # a
    [0.160910839, 0.718971667, 12.9281016],     # b
    [0.302885602, 0.278370726, 0.0152124129],   # c
    [0.170134854, 1.49862703, 0.0283510822],    # d
]

_original_get_element_symbol = kirkland.get_element_symbol


def get_element_symbol(Z):
    if Z == 78:
        return "Pt"
    return _original_get_element_symbol(Z)


kirkland.get_element_symbol = get_element_symbol


# ============================================================
# Load ASE model
# ============================================================

print("=" * 72)
print("QuScope 3-slice transverse-sampling convergence test")
print("=" * 72)

print(f"Model: {MODEL_FILE}")

atom = read(MODEL_FILE)

atomic_numbers = atom.get_atomic_numbers()
symbols = atom.get_chemical_symbols()
positions = atom.get_positions().copy()

composition = Counter(symbols)

print(f"Number of atoms: {len(atom)}")
print(f"Composition: {dict(composition)}")
print(f"Atomic numbers: {sorted(set(atomic_numbers))}")

print("\nASE cell:")
print(atom.cell)

print("\nOriginal atomic extent:")
original_extent = positions.max(axis=0) - positions.min(axis=0)
print(f"x = {original_extent[0]:.3f} Å")
print(f"y = {original_extent[1]:.3f} Å")
print(f"z = {original_extent[2]:.3f} Å")


# ============================================================
# Center particle in a 128 Å transverse FOV
# ============================================================

mins = positions.min(axis=0)
maxs = positions.max(axis=0)
extent = maxs - mins

positions -= mins

positions[:, 0] += (FIELD_OF_VIEW - extent[0]) / 2
positions[:, 1] += (FIELD_OF_VIEW - extent[1]) / 2

print("\nQuScope transverse field:")
print(f"FOV = {FIELD_OF_VIEW:.2f} Å")

print("\nFinal particle bounds:")
print(
    f"x = {positions[:, 0].min():.3f} -> "
    f"{positions[:, 0].max():.3f} Å"
)
print(
    f"y = {positions[:, 1].min():.3f} -> "
    f"{positions[:, 1].max():.3f} Å"
)
print(
    f"z = {positions[:, 2].min():.3f} -> "
    f"{positions[:, 2].max():.3f} Å"
)

print("\nVacuum margins:")
print(
    f"x: {positions[:, 0].min():.3f} Å / "
    f"{FIELD_OF_VIEW - positions[:, 0].max():.3f} Å"
)
print(
    f"y: {positions[:, 1].min():.3f} Å / "
    f"{FIELD_OF_VIEW - positions[:, 1].max():.3f} Å"
)


# ============================================================
# Define physical z slices
# ============================================================

z_min = positions[:, 2].min()
z_max = positions[:, 2].max()

n_slices = int(np.ceil((z_max - z_min) / SLICE_THICKNESS))

slice_indices = np.floor(
    (positions[:, 2] - z_min) / SLICE_THICKNESS
).astype(int)

slice_indices = np.clip(
    slice_indices,
    0,
    n_slices - 1,
)

middle = n_slices // 2

test_slice_indices = np.array(
    [middle - 1, middle, middle + 1],
    dtype=int,
)

print("\nZ slicing:")
print(f"Specimen thickness = {z_max - z_min:.3f} Å")
print(f"Slice thickness    = {SLICE_THICKNESS:.6f} Å")
print(f"Total slices       = {n_slices}")
print(
    "Testing middle slices:",
    test_slice_indices.tolist(),
)


# ============================================================
# Build one potential for one resolution
# ============================================================

def build_test_potential(N, pixel_size):
    """
    Build the same three physical middle slices at a given
    transverse sampling.
    """

    print("\n" + "=" * 72)
    print(f"BUILDING POTENTIAL: {N} x {N}")
    print("=" * 72)

    print(f"Pixel size = {pixel_size:.6f} Å")
    print(f"FOV        = {N * pixel_size:.6f} Å")
    print(f"Slices     = {test_slice_indices.tolist()}")

    xs = np.linspace(
        0.0,
        N * pixel_size,
        N,
        endpoint=False,
    )

    X, Y = np.meshgrid(
        xs,
        xs,
        indexing="ij",
    )

    potentials = []

    for slice_number in test_slice_indices:

        indices = np.where(
            slice_indices == slice_number
        )[0]

        V_slice = np.zeros(
            (N, N),
            dtype=float,
        )

        slice_composition = Counter(
            symbols[index] for index in indices
        )

        z0 = z_min + slice_number * SLICE_THICKNESS
        z1 = z0 + SLICE_THICKNESS

        print(
            f"\nSlice {slice_number}: "
            f"z = {z0:.3f} -> {z1:.3f} Å"
        )
        print(
            f"  atoms       = {len(indices)}"
        )
        print(
            f"  composition = {dict(slice_composition)}"
        )

        for index in indices:

            x_atom = positions[index, 0]
            y_atom = positions[index, 1]
            Z = int(atomic_numbers[index])

            V_atom = kirkland.calculate_2d(
                X,
                Y,
                atom_x=x_atom,
                atom_y=y_atom,
                Z=Z,
            )

            # calculate_2d gives the projected atomic potential
            # used by QuScope's phase-grating convention.
            V_slice += V_atom

        potentials.append(V_slice)

        print(
            f"  potential min = {V_slice.min():.6e}"
        )
        print(
            f"  potential max = {V_slice.max():.6e}"
        )

    return np.asarray(potentials)


# ============================================================
# Run QuScope for each resolution
# ============================================================

for N, pixel_size in RESOLUTIONS.items():

    result_file = (
        OUTPUT_DIR
        / f"Pt5nm_middle3_{N}x{N}.npz"
    )

    potential_file = (
        OUTPUT_DIR
        / f"Pt5nm_middle3_potential_{N}x{N}.npz"
    )

    # --------------------------------------------------------
    # Skip an already completed resolution
    # --------------------------------------------------------

    if result_file.exists():

        print("\n" + "=" * 72)
        print(f"RESULT ALREADY EXISTS: {result_file}")
        print("Skipping this resolution.")
        print("=" * 72)

        continue

    # --------------------------------------------------------
    # Build potential
    # --------------------------------------------------------

    t0 = time.perf_counter()

    potentials = build_test_potential(
        N,
        pixel_size,
    )

    potential_time = time.perf_counter() - t0

    print(
        f"\nPotential generation time: "
        f"{potential_time / 60:.2f} min"
    )

    # Save potential before the expensive quantum calculation.
    np.savez_compressed(
        potential_file,
        potentials=potentials,
        grid_size=N,
        pixel_size=pixel_size,
        field_of_view=FIELD_OF_VIEW,
        slice_thickness=SLICE_THICKNESS,
        slice_indices=test_slice_indices,
        acceleration_voltage=ACCELERATION_VOLTAGE,
        model_file=str(MODEL_FILE),
    )

    print(f"Saved potential to: {potential_file}")

    # --------------------------------------------------------
    # Build QuScope circuit
    # --------------------------------------------------------

    params = QuantumMultisliceParameters(
        acceleration_voltage=ACCELERATION_VOLTAGE,
        grid_size=N,
        pixel_size=pixel_size,
        defocus=DEFOCUS,
        cs=CS,
        slice_thickness=SLICE_THICKNESS,
    )

    print("\nBuilding QuScope circuit...")

    t0 = time.perf_counter()

    circuit = QuantumMultisliceCircuit(params)

    circuit_time = time.perf_counter() - t0

    print(
        f"Circuit construction: "
        f"{circuit_time:.2f} s"
    )

    # --------------------------------------------------------
    # Run quantum multislice
    # --------------------------------------------------------

    print("\n" + "=" * 72)
    print(
        f"STARTING {N} x {N} "
        f"QUANTUM MULTISLICE"
    )
    print("=" * 72)

    t0 = time.perf_counter()

    result = circuit.simulate(
        potentials
    )

    simulation_time = time.perf_counter() - t0

    print(
        f"\nSimulation time: "
        f"{simulation_time / 60:.2f} min"
    )

    # --------------------------------------------------------
    # Save result immediately
    # --------------------------------------------------------

    np.savez_compressed(
        result_file,
        statevector=result["statevector"],
        wave_function=result["wave_function"],
        amplitude=result["amplitude"],
        phase=result["phase"],
        intensity=result["intensity"],
        grid_size=N,
        pixel_size=pixel_size,
        field_of_view=FIELD_OF_VIEW,
        slice_thickness=SLICE_THICKNESS,
        slice_indices=test_slice_indices,
        acceleration_voltage=ACCELERATION_VOLTAGE,
        defocus=DEFOCUS,
        cs=CS,
        model_file=str(MODEL_FILE),
        potential_file=str(potential_file),
    )

    print(f"Saved result to: {result_file}")

    print("\nResult:")
    print(
        f"  wave function shape = "
        f"{result['wave_function'].shape}"
    )
    print(
        f"  amplitude max       = "
        f"{result['amplitude'].max():.6e}"
    )
    print(
        f"  intensity sum       = "
        f"{result['intensity'].sum():.6f}"
    )

    print(
        f"\nCompleted {N} x {N} resolution."
    )


print("\n" + "=" * 72)
print("ALL REQUESTED CONVERGENCE CALCULATIONS COMPLETE")
print("=" * 72)
