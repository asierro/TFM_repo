import sisl
import numpy as np
import matplotlib.pyplot as plt

zgnr8 = sisl.geom.zgnr(4).tile(38, 0).move([0, -10., 0])

sc_x = zgnr8.cell[0, 0]
sc_y = 8.52
sc_z = 50.
zgnr8.set_supercell([sc_x, sc_y, sc_z])

dp_pp = zgnr8.remove([3,
                      48, 52, 61, 65, 74, 78,
                      127, 131,
                      136, 140, 149, 153,
                      277, 281, 290, 294, 303,
                      224, 228,
                      202, 206, 215, 219])

H = sisl.Atom(1)
H_list = [0, 3, 7,
          6, 10,
          42, 46, 49,
          39, 43, 50, 54, 61, 65, 72, 76, 80, 79, 75, 68, 64, 53, 57,
          116, 112, 111, 115, 119, 126, 129,
          120, 123, 130, 134, 141, 145, 149, 148, 144, 137, 133,
          196, 192, 185, 181, 180, 184, 188, 195, 199, 206, 209, 200, 203, 210, 214, 218, 217, 213,
          276, 272, 265, 261, 254, 250, 249, 253, 257, 264, 268, 275, 279]
dp_pp.atoms.replace(H_list, H)

#dp_pp = dp_pp.tile(2, 1).tile(2,0)

# dp_pp.write('STRUCT_pp_ortho.fdf')

sisl.plot(dp_pp, s=5); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
plt.tight_layout(); plt.show()

# sisl.plot(dp_pp, s=5); ax = plt.gca(); ax.set_aspect('equal', adjustable='box'); plt.tight_layout()
# fig=plt.gcf(); w, h = fig.get_size_inches(); fig.set_size_inches(w * 2, h * 2); plt.savefig('final.pdf')
