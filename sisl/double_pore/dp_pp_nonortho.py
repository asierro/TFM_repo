import sisl
import numpy as np
import matplotlib.pyplot as plt

zgnr8 = sisl.geom.zgnr(8).tile(10, 0).rotate(180, [0, 0, 1])
zgnr8 = zgnr8.move([-zgnr8.axyz()[159][0], 24.91, 0])

coords = zgnr8.axyz()

a = 1.42
sc_1 = [19 / 2 * np.sqrt(3) * a, 9 / 2 * a, 0]
sc_2 = [0., 6 * a, 0]
sc_3 = [0., 0., 50.]
zgnr8.set_supercell([sc_1, sc_2, sc_3])

remove_list = []
for ia in range(zgnr8.na):
    xa = coords[ia][0]
    ya = coords[ia][1]
    if  ya < 9 / (19 * np.sqrt(3)) * xa:
        remove_list.append(ia)
    elif ya > 6 * a + 9 / (19 * np.sqrt(3)) * xa:
        remove_list.append(ia)

zgnr8 = zgnr8.remove(remove_list)

dp_pp = zgnr8.remove([0, 1, 2, 3, 4,
                      7, 15, 19, 28, 31,
                      79])

H = sisl.Atom(1)
H_list = [68, 64, 60, 61, 65,
          23,
          4, 8, 15, 19, 26, 29, 30, 22, 16, 12, 5, 2]
dp_pp.atoms.replace(H_list, H)

# dp_pp.write('STRUCT_pp.fdf')

# plt.gca().set_aspect('equal', adjustable='box'); sisl.plot(dp_pp); plt.show()

# dp_pp = dp_pp.tile(2, 1).tile(2,0)
# plt.gca().set_aspect('equal', adjustable='box'); sisl.plot(dp_pp, s=20); plt.show()
