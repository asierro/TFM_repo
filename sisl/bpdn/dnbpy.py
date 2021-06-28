import sisl
import numpy as np
import matplotlib.pyplot as plt

'''
Creation of dinitro bypiridine
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


ang_list1 = [116.6, -117.8, -123.8]
dist_list1 = [1.337, 1.333, 1.383]
pos_list = pos_gen(ang_list1, dist_list1)

ang_list2 = [-120.6, 121.7]
dist_list2 = [1.372, 1.368]

pos_list.extend(pos_gen(ang_list2, dist_list2))

pos_list.append(np.zeros(3))

# H atoms
pos_list.append(next_coord(pos_list[1], -61.2, 1.103))
pos_list.append(next_coord(pos_list[4], 61.1, 1.103))

# ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# plt.scatter([pos[0] for pos in pos_list], [pos[1] for pos in pos_list]); plt.show()

H = sisl.Atom(1)
C = sisl.Atom(6)
N = sisl.Atom(7)
O = sisl.Atom(8)

# Geometry object for the pyridil (w/o nitro group)
a_list=[N, C, C, C, C, C, H, H]
pyridyl_atoms = sisl.Atoms(atoms=a_list, na=8)
pyridyl_sc = sisl.SuperCell([50., 50., 50.,]) # TODO: set to this at end
# pyridyl_sc = sisl.SuperCell([15., 6., 50.,])
pyridyl = sisl.Geometry(pos_list, atoms=pyridyl_atoms, sc=pyridyl_sc)

# Geometry object for the nitro grup
a_list2=[N, O, O]
nitro_atoms = sisl.Atoms(atoms=a_list2, na=3)
nitro_sc = sisl.SuperCell([5., 5., 5.,])
x = 1.219 * np.cos(27.8 * np.pi / 180)
y = 1.219 * np.sin(27.8 * np.pi / 180)
nitro_coords = [np.array([0. ,0. ,0]),
                np.array([-x, y, 0.]),
                np.array([x, y, 0.])]
nitro = sisl.Geometry(nitro_coords, atoms=nitro_atoms, sc=nitro_sc)
nitro = nitro.rotate(30, [0., 1., 0.])
nitro = nitro.rotate(29.9, [0., 0., 1.])

# ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# sisl.plot(nitro); plt.show()

offset = next_coord(pyridyl.axyz()[3], 119.9, 1.479)
nitro_pyridyl = pyridyl.add(nitro, offset=offset)

# ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# sisl.plot(nitro_pyridyl); plt.show()

# Geometry object for the remaining part of the oligophenylene-ethynylene
ang_list3 = [0., 0., 0., 120., -120., -120., -120., -120.]
dist_list3 = [1.45, 1.18, 1.45, 1.41, 1.41, 1.41, 1.41, 1.41]

ope_coords = pos_gen(ang_list3, dist_list3)
# H atoms
ope_coords.append(next_coord(ope_coords[3], -120., 1.103))
ope_coords.append(next_coord(ope_coords[4], -60., 1.103))
ope_coords.append(next_coord(ope_coords[6], 60., 1.103))
ope_coords.append(next_coord(ope_coords[7], 120., 1.103))
# HERE I DEFINE THE END BIT
# in this case an H atoms
# Remember to change a_list3 and na accordingly
ope_coords.append(next_coord(ope_coords[5], 0., 1.103))

#
# ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# plt.scatter([pos[0] for pos in ope_coords], [pos[1] for pos in ope_coords]); plt.show()

a_list3 = [C, C, C, C, C, C, C, C, H, H, H, H, H]
ope_atoms = sisl.Atoms(atoms=a_list3, na=13)
ope_sc = sisl.SuperCell([5., 5., 5.,])
nitro = sisl.Geometry(ope_coords, atoms=ope_atoms, sc=ope_sc)

half_bpdn = nitro_pyridyl.add(nitro, offset=nitro_pyridyl.axyz()[2])

# ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# sisl.plot(half_bpdn); plt.show()

# Make other half
half2_bpdn = half_bpdn.rotate(180, [0., 0., 1.,])
half2_bpdn = half2_bpdn.rotate(-1.08, [1., 0., 0.,], rad=True)
half2_bpdn = half2_bpdn.rotate(-30, [0., 1., 0.,])


# aan = half2_bpdn.axyz()

# FINAL GEOMETRY

z = -1.492 * np.sin(30 * np.pi / 180)
x = -1.492 * np.cos(30 * np.pi / 180)
bpdn = half_bpdn.add(half2_bpdn, offset=[x, 0., z])

# ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# sisl.plot(bpdn); plt.show()

sisl.get_sile('test_bpdn_ope.xsf', 'w').write_geometry(bpdn)
bpdn.write('STRUCT_bpdn_ope.fdf')
