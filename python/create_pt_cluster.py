from pathlib import Path
import logging
import sys

from ase.cluster import wulff_construction
from ase.io import write

ROOT = Path("/zhome/f1/0/167759/10387")
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pt_cluster.functions import Zoneaxis4, stable_solution

DATA = ROOT / "data"
MODEL_DIR = DATA / "models"
LOG_DIR = DATA / "logs"

for folder in [MODEL_DIR, LOG_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    filename=LOG_DIR / "create_pt_cluster.log",
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

surfaces = [
    (1, 0, 0), (1, 1, 0), (1, 1, 1),
    (3, 3, 2), (3, 2, 2), (2, 2, 1),
    (3, 3, 1), (2, 1, 1), (3, 2, 1),
    (3, 1, 1), (3, 1, 0), (2, 1, 0), (3, 2, 0),
]

esurf = [
    0.116, 0.117, 0.093,
    0.097, 0.099, 0.100,
    0.106, 0.110, 0.110,
    0.112, 0.117, 0.118, 0.118,
]

lc = 3.9236
size = 4045
zone_axis = (1, 1, 0)

Pt = wulff_construction(
    "Pt",
    surfaces,
    esurf,
    size=size,
    structure="fcc",
    rounding="closest",
    latticeconstant=lc,
)

buffer = Pt.get_diameter(method="shape")
print(f"diameter = {buffer}", f"atoms = {len(Pt)}")
logging.info("Raw Pt cluster: diameter=%s A, atoms=%s", buffer, len(Pt))

write(MODEL_DIR / "Pt5nm_prime.traj", Pt)
write(MODEL_DIR / "Pt5nm_prime.xyz", Pt)

atom = Pt.copy()
atom = stable_solution(atom)
theta, phi = Zoneaxis4(zone_axis)

atom.rotate(theta, (1, 0, 0), center=(0, 0, 0), rotate_cell=False)
atom.rotate(phi, (0, 1, 0), center=(0, 0, 0), rotate_cell=False)

atom.center(buffer)

final_model = MODEL_DIR / "Pt5nm_cov.traj"
final_model_xyz = MODEL_DIR / "Pt5nm_cov.xyz"

write(final_model, atom)
write(final_model_xyz, atom)

logging.info("Zone axis %s rotation: theta=%s, phi=%s", zone_axis, theta, phi)
logging.info("Final cell lengths: %s", atom.cell.lengths())
logging.info("Saved final model to %s", final_model)

print(f"saved final model: {final_model}")
print(f"cell lengths = {atom.cell.lengths()}")