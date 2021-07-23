import sisl
import numpy as np
import matplotlib.pyplot as plt

zgnr8 = sisl.geom.zgnr(4, bond=1.44).tile(19, 0).move([0, -10., 0])

sc_x = zgnr8.cell[0, 0]
sc_y = 8.52
sc_z = 50.
zgnr8.set_supercell([sc_x, sc_y, sc_z])

dp_pp = zgnr8.remove([48, 52, 61, 65, 74, 78,
                      126, 130, 137, 141, 148, 0])

H = sisl.Atom(1)
H_list = [2, 6, 10, 7, 3,
          39, 43, 42, 46, 49, 50, 54, 61, 65, 72, 76, 80, 79, 75, 68, 64, 57, 53,
          133, 130, 124, 121, 115, 111, 112, 116, 119, 125, 128, 134, 137]
dp_pp.atoms.replace(H_list, H)

dp_pp = dp_pp.tile(2, 1).tile(2,0)

# dp_pp.write('STRUCT_pp_alter.fdf')

coords = dp_pp.axyz()
ats = list(dp_pp.iter_species())

for ia, a, idx_specie in dp_pp.iter_species():
    if idx_specie == 1:  # H to move
        close = dp_pp.close(ia, R=1.45)
        for ja in close: # C to move towards
            if ats[ja][2] == 0:
                dp_pp = dp_pp.move((1.44 - 1.105) * (coords[ja] - coords[ia]) / 1.44, atoms=ia)

# by using ax.set_xlim, y_lim and atom_indices = True determine the position of the atom to move to (dp_pp.move),
# then redefine the sc (calculate size mathematically) and remove outliers with the function I wrote.

f = plt.figure()
f.clear()
plt.close(f)
plt.clf()

sisl.plot(dp_pp, s=10); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
plt.tight_layout(); plt.show()

# sisl.plot(dp_pp, s=5, atom_indices=True); ax = plt.gca(); ax.set_aspect('equal', adjustable='box'); plt.tight_layout()
# fig=plt.gcf(); w, h = fig.get_size_inches(); fig.set_size_inches(w * 2, h * 2); plt.savefig('alter_h.pdf')
