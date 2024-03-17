from inertia_moment import Inertia
from os.path import join

BASE_PATH = r"E:\C_backup\Documents-2\Documents\Simulations\UBI\room_temp"
SYSTEM_NAME = "60mM"

top = join(BASE_PATH, SYSTEM_NAME, "centered.gro")
trj = join(BASE_PATH, SYSTEM_NAME, "centered.xtc")

m = Inertia(top, trj)

m.select_atoms("(resname SDS and name OS1) or protein")
# m.run()

m.run(start=1800, end=2350, skip=1)
shape = m.determine_shape()

print(Inertia.shape_mapping[shape])
