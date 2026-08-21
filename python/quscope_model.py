import numpy as np
from ase.io import read
from collections import Counter
from quscope.ctem import KirklandPotential

# ============================================================
# Load Pt model
# ============================================================

atom = read(data / "models/Pt5nm_cov.traj")

positions = atom.get_positions().copy()


# ============================================================
# QuScope grid
# ============================================================

N = 256
pixel_size = 0.25  # Å/pixel

xs = np.linspace(
    0,
    pixel_size * N,
    N,
    endpoint=False,
)

X, Y = np.meshgrid(
    xs,
    xs,
    indexing="ij",
)

field_of_view = N * pixel_size

print(f"Field of view: {field_of_view:.2f} Å")


# ============================================================
# Center particle in QuScope transverse field of view
# ============================================================

positions = atom.get_positions().copy()

# Original particle bounds
mins = positions.min(axis=0)
maxs = positions.max(axis=0)

extent = maxs - mins

print("\nOriginal particle extent:")
print(f"x = {extent[0]:.2f} Å")
print(f"y = {extent[1]:.2f} Å")
print(f"z = {extent[2]:.2f} Å")


# ------------------------------------------------------------
# Shift particle so its minimum x/y/z coordinates are zero
# ------------------------------------------------------------

positions -= mins


# ------------------------------------------------------------
# Center x/y in the QuScope field
# ------------------------------------------------------------

positions[:, 0] += (
    field_of_view - extent[0]
) / 2

positions[:, 1] += (
    field_of_view - extent[1]
) / 2


# z remains:
# 0 → specimen thickness


# ============================================================
# Check final coordinates
# ============================================================

print("\nFinal particle coordinates:")

print(
    f"x = {positions[:, 0].min():.2f}"
    f" → {positions[:, 0].max():.2f} Å"
)

print(
    f"y = {positions[:, 1].min():.2f}"
    f" → {positions[:, 1].max():.2f} Å"
)

print(
    f"z = {positions[:, 2].min():.2f}"
    f" → {positions[:, 2].max():.2f} Å"
)

print("\nFinal extent:")
final_extent = positions.max(axis=0) - positions.min(axis=0)

print(f"x = {final_extent[0]:.2f} Å")
print(f"y = {final_extent[1]:.2f} Å")
print(f"z = {final_extent[2]:.2f} Å")


# ============================================================
# Kirkland Pt potential
# ============================================================

kirkland = KirklandPotential()

# Pt Kirkland parameterization used in the QuScope calculation.
# Source: E. J. Kirkland, Advanced Computing in Electron Microscopy,
# 2nd ed., Springer (2010), as provided through the abTEM parameter data.
kirkland.params_dict["Pt"] = [
    [0.98469794, 2.73987079, 3.61696715],
    [0.160910839, 0.718971667, 12.9281016],
    [0.302885602, 0.278370726, 0.0152124129],
    [0.170134854, 1.49862703, 0.0283510822],
]

original_get_element_symbol = kirkland.get_element_symbol

def get_element_symbol(Z):
    if Z == 78:
        return "Pt"
    return original_get_element_symbol(Z)

kirkland.get_element_symbol = get_element_symbol


# ============================================================
# Slice along z
# ============================================================

slice_thickness = 1.0  # Å

z_min = positions[:, 2].min()
z_max = positions[:, 2].max()

n_slices = int(
    np.ceil((z_max - z_min) / slice_thickness)
)

print(f"\nSpecimen thickness: {z_max - z_min:.2f} Å")
print(f"Slice thickness:    {slice_thickness:.2f} Å")
print(f"Number of slices:   {n_slices}")


# ============================================================
# Assign atoms to slices
# ============================================================

slice_indices = np.floor(
    (positions[:, 2] - z_min) / slice_thickness
).astype(int)

slice_indices = np.clip(
    slice_indices,
    0,
    n_slices - 1,
)


# ============================================================
# Build potentials
# ============================================================

potentials = []
atomic_numbers = atom.get_atomic_numbers()

from collections import Counter

for i in range(n_slices):

    indices = np.where(
        slice_indices == i
    )[0]

    V_slice = np.zeros(
        (N, N),
        dtype=float,
    )

    for index in indices:

        x_atom = positions[index, 0]
        y_atom = positions[index, 1]
        Z = atomic_numbers[index]

        V_atom = kirkland.calculate_2d(
            X,
            Y,
            atom_x=x_atom,
            atom_y=y_atom,
            Z=Z,
        )

        V_slice += V_atom

    potentials.append(V_slice)

    # Composition of THIS slice
    slice_composition = Counter(
        symbols[index] for index in indices
    )

    print(
        f"Slice {i+1:3d}/{n_slices}: "
        f"{len(indices):3d} atoms | "
        f"{dict(slice_composition)} | "
        f"max = {V_slice.max():.4e}"
    )

potentials = np.asarray(potentials)

print("\nFinal potential:")
print("shape =", potentials.shape)