import sisl
import numpy as np
import matplotlib.pyplot as plt

'''
Creation of 2-2'-dinitro biphenyl dithiol.
(Based on bpdn-ope.py)
'''


def next_coord(pos_prev, angle, dist):
    # angle: wrt +x axis
    pos = np.zeros(3)
    pos[0] = pos_prev[0] + dist * np.cos(angle * np.pi / 180)
    pos[1] = pos_prev[1] + dist * np.sin(angle * np.pi / 180)
    return pos


def pos_gen(ang_list, dist_list):
    # ang: wrt last bond
    assert len(ang_list) == len(dist_list)
    pos_list = []
    angle = 0.
    for i, ang in enumerate(ang_list):
        angle += ang - np.sign(ang) * 180
        if i != 0:
            pos_list.append(next_coord(pos_list[-1], angle, dist_list[i]))
        else:
            pos_list.append(next_coord(np.zeros(3), angle, dist_list[i]))
    return pos_list


H = sisl.Atom(1)
C = sisl.Atom(6)
N = sisl.Atom(7)
O = sisl.Atom(8)
S = sisl.Atom(16)

# Geometry object for the phenyl-thiol
ang_list = [0., 120., -120., -120., -120., -120.]
dist_list = [0., 1.41, 1.41, 1.41, 1.41, 1.41]

phenyl_coords = pos_gen(ang_list, dist_list)
# H atoms
phenyl_coords.append(next_coord(phenyl_coords[1], -120., 1.103))
phenyl_coords.append(next_coord(phenyl_coords[2], -60., 1.103))
phenyl_coords.append(next_coord(phenyl_coords[4], 60., 1.103))
# HERE I DEFINE THE END BIT
# in this case an H atom
# Remember to change a_list and na accordingly
phenyl_coords.append(next_coord(phenyl_coords[3], 0., 1.8))
phenyl_coords.append(next_coord(phenyl_coords[9], 0., 1.34))

#
# ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# plt.scatter([pos[0] for pos in ope_coords], [pos[1] for pos in ope_coords]); plt.show()

a_list = [C, C, C, C, C, C, H, H, H, S, H]
phenyl_atoms = sisl.Atoms(atoms=a_list, na=11)
phenyl_sc = sisl.SuperCell([20., 10., 50.])
phenyl = sisl.Geometry(phenyl_coords, atoms=phenyl_atoms, sc=phenyl_sc)

# ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# sisl.plot(phenyl); plt.show()

# Geometry object for the nitro grup
a_list2=[N, O, O]
nitro_atoms = sisl.Atoms(atoms=a_list2, na=3)
nitro_sc = sisl.SuperCell([5., 5., 5.,])
x = 1.219 * np.cos(27.8 * np.pi / 180)
y = 1.219 * np.sin(27.8 * np.pi / 180)
nitro_coords = [np.array([0., 0., 0]),
                np.array([-x, y, 0.]),
                np.array([x, y, 0.])]
nitro = sisl.Geometry(nitro_coords, atoms=nitro_atoms, sc=nitro_sc)
# nitro = nitro.rotate(30, [0., 1., 0.])
nitro = nitro.rotate(30., [0., 0., 1.])

# ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# sisl.plot(nitro); plt.show()

# Appending the nitro group
offset = next_coord(phenyl.axyz()[5], 119.9, 1.479)
nitro_phenyl = phenyl.add(nitro, offset=offset)

# ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# sisl.plot(nitro_phenyl); plt.show()

# 2nd half
half2 = nitro_phenyl.rotate(180, [0., 0., 1.])
half2 = half2.rotate(90, [1., 0., 0.])

dbdt = nitro_phenyl.add(half2, offset=[-1.492, 0., 0.])
dbdt = dbdt.move([10., 3., 3.])

# ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# sisl.plot(dbdt); plt.show()
#
# sisl.get_sile('test_dbdt.xsf', 'w').write_geometry(dbdt)
# dbdt.write('STRUCT_dbdt.fdf')
