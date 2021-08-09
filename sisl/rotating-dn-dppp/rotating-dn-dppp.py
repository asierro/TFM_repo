import sisl
import numpy as np
import matplotlib.pyplot as plt

'''
Creation of dp_pp structure functionalized with two nitro groups, replacing the biphenyl bridges with
2,2'-dinitrobiphenyl, with different twist angles of the phenyl groups (while keeping nitro groups planar wrt the 
corresponding ring.
Based off of the relaxed planar structure.
'''


def next_coord(pos_prev, angle, dist):
    # angle: wrt +x axis
    pos = np.zeros(3)
    pos[0] = pos_prev[0] + dist * np.cos(angle * np.pi / 180)
    pos[1] = pos_prev[1] + dist * np.sin(angle * np.pi / 180)
    pos[2] = pos_prev[2]
    return pos


XV_sile = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/dp_pp/relax/NPG_PP-nonplanar.XV')
dp_pp_ortho = sisl.Geometry.read(XV_sile)

coords = dp_pp_ortho.axyz()

# Move H atoms belonging to phenyls that were on the other side of the supercell
sc_y = dp_pp_ortho.cell[1]
dp_pp_ortho = dp_pp_ortho.move(-sc_y, atoms=19).move(-sc_y, atoms=24)

phenyl1 = [103, 104, 105, 106,
           16, 18, 19, 13]
phenyl2 = [109, 110, 111, 112,
           0, 17, 14, 10]
phenyl3 = [119, 118, 122, 116,
           22, 23, 26, 21]
phenyl4 = [114, 121, 91, 125,
           25, 34, 27, 24]

planar_dppp = dp_pp_ortho.rotate(-28.67, coords[102] - coords[107], atoms=phenyl1, origo=coords[107], only='xyz') \
                         .rotate(28.23, coords[108] - coords[113], atoms=phenyl2, origo=coords[113], only='xyz') \
                         .rotate(-28.28, coords[115] - coords[117], atoms=phenyl3, origo=coords[117], only='xyz') \
                         .rotate(28.85, coords[90] - coords[120], atoms=phenyl4, origo=coords[120], only='xyz')

# sisl.plot(planar_dppp, s=10, atom_indices=True); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# plt.tight_layout(); plt.show()

# Geometry object for the nitro group
N = sisl.Atom(7)
O = sisl.Atom(8)
a_list = [N, O, O]
nitro_atoms = sisl.Atoms(atoms=a_list, na=3)
nitro_sc = sisl.SuperCell([5., 5., 5.,])
x = 1.219 * np.cos(27.8 * np.pi / 180)
y = 1.219 * np.sin(27.8 * np.pi / 180)
nitro_coords = [np.array([0., 0., 0]),
                np.array([-x, y, 0.]),
                np.array([x, y, 0.])]
nitro = sisl.Geometry(nitro_coords, atoms=nitro_atoms, sc=nitro_sc)
ud_nitro = nitro.rotate(180., [0., 0., 1.])  # Upside-down nitro

# ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# sisl.plot(ud_nitro); plt.show()

# Remove H atoms where nitro groups will go
planar_dppp = planar_dppp.remove([18, 14, 21, 25])

# Add nitro groups
offset1 = next_coord(planar_dppp.axyz()[100], 90, 1.483)
offset2 = next_coord(planar_dppp.axyz()[107], -90, 1.483)
offset3 = next_coord(planar_dppp.axyz()[112], -90, 1.483)
offset4 = next_coord(planar_dppp.axyz()[110], 90, 1.483)
dn_dppp = planar_dppp.add(nitro, offset=offset1).add(ud_nitro, offset=offset2).add(ud_nitro, offset=offset3) \
                     .add(nitro, offset=offset4)

coords = dn_dppp.axyz()

nphenyl1 = [99, 100, 101, 102,  # C
           15, 17, 13,         # H
           136, 137, 138]      # Nitro
nphenyl2 = [105, 106, 107, 108,
           10, 16, 0,
           139, 140, 141]
nphenyl3 = [112, 114, 115, 118,
           22, 19, 20,
           142, 143, 144]
nphenyl4 = [87, 121, 117, 110,
           23, 21, 30,
           145, 146, 147]

# dn_dppp = dn_dppp.rotate(90, coords[103] - coords[98], atoms=nphenyl1, origo=coords[98], only='xyz')
# dn_dppp = dn_dppp.rotate(90, coords[109] - coords[104], atoms=nphenyl2, origo=coords[104], only='xyz')
# dn_dppp = dn_dppp.rotate(90, coords[113] - coords[111], atoms=nphenyl3, origo=coords[111], only='xyz')
# dn_dppp = dn_dppp.rotate(90, coords[116] - coords[86], atoms=nphenyl4, origo=coords[86], only='xyz')

# sisl.plot(dn_dppp, s=10); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# plt.tight_layout(); plt.show()


def rotate_nphenyl(geom, which, angle):
    '''
    Rotate a phenyl from the planar position, the rotation vectors go from right to left.
    geom must be based on dp_pp_ortho (same indices)
    '''
    assert which in ['1', '2', '3', '4']
    if which == '1':
        rotated = geom.rotate(angle, coords[103] - coords[98], atoms=nphenyl1, origo=coords[98], only='xyz')
    elif which == '2':
        rotated = geom.rotate(angle, coords[109] - coords[104], atoms=nphenyl2, origo=coords[104], only='xyz')
    elif which == '3':
        rotated = geom.rotate(angle, coords[113] - coords[111], atoms=nphenyl3, origo=coords[111], only='xyz')
    else:
        rotated = geom.rotate(angle, coords[116] - coords[86], atoms=nphenyl4, origo=coords[86], only='xyz')
    return rotated


def rotate_dn_dppp(geom, ang1, ang2, ang3, ang4):
    rp = rotate_nphenyl(geom, '1', ang1)
    rp = rotate_nphenyl(rp, '2', ang2)
    rp = rotate_nphenyl(rp, '3', ang3)
    return rotate_nphenyl(rp, '4', ang4)


ang_a = 30
for ang_a in range(0, 100, 10):
    for ang_b in range(0, 180, 10):
        print(ang_a, ang_b)
        rotating_dn_dppp = rotate_dn_dppp(dn_dppp, ang_a, ang_a + ang_b, -ang_a, -ang_a - ang_b)
        # rotating_dppp.write('structures/STRUCT_{}_{}.fdf'.format(ang_a, ang_b))
        f = plt.figure()
        f.clear()
        plt.close(f)
        plt.clf()
        sisl.plot(rotating_dn_dppp, s=10); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
        plt.tight_layout()
        plt.savefig(
            '/home/asier/Documents/master/tfm/sisl/rotating-dn-dppp/gif/rot_{0:03d}_{1:03d}.png'.format(ang_a, ang_b))
