import sisl
import numpy as np
import matplotlib.pyplot as plt
'''
Prepare an orthogonal DP PP for relaxation, based on the relaxed NPG
'''


fdf = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/npg/4_bands/RUN.fdf')

npg = sisl.Geometry.read(fdf, output=True)  # read from output files

sc_y = npg.cell[1, 1]

npg = npg.tile(2, 1)

a = 1.43
sc_x = npg.cell[0, 0] + 12 * a * np.sqrt(3) / 2
sc_z = 50.
npg.set_supercell([sc_x, sc_y, sc_z])


pos = npg.axyz()

# ---------------------------------------------------

r_atoms = []
l_atoms = []
for ia, at, idx_specie in npg.iter_species():
    if idx_specie == 0:  # carbon atoms
        if pos[ia][0] > 25.:
            r_atoms.append(ia)
        if pos[ia][0] < 7.7:
            l_atoms.append(ia)

r_h = []
for ra in r_atoms:
    idx = npg.close(ra, [0.1, 1.2])  # find hydrogen atoms
    if idx[1].size > 0:
        r_h.append(idx[1][0])
r_atoms.extend(r_h)

l_h = []
for la in l_atoms:
    idx = npg.close(la, [0.1, 1.2])  # find hydrogen atoms
    if idx[1].size > 0:
        l_h.append(idx[1][0])
l_atoms.extend(l_h)

c_atoms = list(set(npg.iter()).difference(set(r_atoms), set(l_atoms)))

# ---------------------------------------------------

x = 6 * a * np.sqrt(3.) / 2
y = 6 * a / 2

npg = npg.move([x, 0., 0.], atoms=c_atoms)
npg = npg.move([2 * x, 0., 0.], atoms=r_atoms)
npg = npg.move([0., -y, 0.], atoms=c_atoms)

# ---------------------------------------------------

# H = sisl.Atom(1)
# C = sisl.Atom(6)  # change for original Atom objects?

H = npg.atoms.atom[1]
C = npg.atoms.atom[0]

a_list=[C, C, C, C, C, C, H, H, H, H]
phenyl_atoms = sisl.Atoms(atoms=a_list, na=10)
s3 = np.sqrt(3.)
b = 1.1033
phenyl_c = a * np.array([[0.5 * s3, -0.5, 0.],
                         [s3, 0., 0.],
                         [1.5 * s3, -0.5, 0.],
                         [1.5 * s3, -1.5, 0.],
                         [s3, -2., 0.],
                         [0.5 * s3, -1.5, 0.]])

phenyl_h = np.zeros([4,3])
phenyl_h[0] = phenyl_c[1] + b * np.array([0., 1., 0.])
phenyl_h[1] = phenyl_c[2] + b * np.array([0.5 * s3, 0.5, 0.])
phenyl_h[2] = phenyl_c[4] + b * np.array([0., -1., 0.])
phenyl_h[3] = phenyl_c[5] + b * np.array([-0.5 * s3, -0.5, 0.])

phenyl_coords = list(phenyl_c)
phenyl_coords.extend(list(phenyl_h))
phenyl_coords = np.array(phenyl_coords)

phenyl_sc = sisl.SuperCell([4., 4., 50.,])

phenyl = sisl.Geometry(phenyl_coords, atoms=phenyl_atoms, sc=phenyl_sc)

# ------------------------------------------------

npg = npg.add(phenyl, offset=pos[22])

pos = npg.axyz()

npg = npg.add(phenyl, offset=pos[203])

pos = npg.axyz()

# ------------------------------------------------

flipped_coords = phenyl_coords
for i in range(10):
    flipped_coords[i][1] *= -1.
flipped = sisl.Geometry(flipped_coords, atoms=phenyl_atoms, sc=phenyl_sc)

npg = npg.add(flipped, offset=pos[69])

pos = npg.axyz()

npg = npg.add(flipped, offset=pos[223])

# ------------------------------------------------

pos = npg.axyz()

remove_list = []
for ia in npg.iter():
    if pos[ia][1] < 0. or pos[ia][1] > sc_y:
        remove_list.append(ia)

npg = npg.remove(remove_list)

# npg = npg.tile(2, 1).tile(2,0)
# npg.write('STRUCT_pre.fdf')
# ------------------------------------------------

# f = plt.figure()
# f.clear()
# plt.close(f)
# plt.clf()
#
sisl.plot(npg, s=10); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
plt.tight_layout(); plt.show()
