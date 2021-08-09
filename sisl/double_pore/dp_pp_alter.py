import sisl
import numpy as np
import matplotlib.pyplot as plt


def remove_outliers(geom):
    '''
    This routine removes all atoms in a geometry which are located outside the supercell (the parallepiped defined by
    the sc vectors)
    '''
    sc = geom.cell
    axyz = geom.axyz()
    # deal with the signs in a more systematic way (here I had to change them from last time)
    u, v, w = (np.cross(sc[0], sc[1]), -np.cross(sc[0], sc[2]), np.cross(sc[1], sc[2]))
    na = axyz.shape[0]  # better ways
    udotaxyz = np.einsum('ij, ij->i', np.tile(u, (na, 1)), axyz)
    vdotaxyz = np.einsum('ij, ij->i', np.tile(v, (na, 1)), axyz)
    wdotaxyz = np.einsum('ij, ij->i', np.tile(w, (na, 1)), axyz)
    c1 = np.logical_and(0 < udotaxyz, udotaxyz <= np.dot(u, sc[2]))
    c2 = np.logical_and(0 < vdotaxyz, vdotaxyz <= np.dot(v, sc[1]))
    c3 = np.logical_and(0 < wdotaxyz, wdotaxyz <= np.dot(w, sc[0]))
    in_bool = np.logical_not(np.logical_and(c1, np.logical_and(c2, c3)))
    rem_list = [np.where(in_bool)[0]]  # list not necessary?
    return geom.remove(rem_list)


zgnr8 = sisl.geom.zgnr(4, bond=1.44).tile(19, 0).move([0, -10., 0])

sc_x = zgnr8.cell[0, 0]
sc_y = 12 / 2 * 1.44
sc_z = 50.
zgnr8.set_supercell([sc_x, sc_y, sc_z])

dp_pp = zgnr8.remove([48, 52, 61, 65, 74, 78,
                      126, 130, 137, 141, 148, 0])

H = sisl.Atom(1)
H_list = [2, 6, 10, 7, 3,
          39, 43, 42, 46, 49, 50, 54, 61, 65, 72, 76, 80, 79, 75, 68, 64, 57, 53,
          133, 130, 124, 121, 115, 111, 112, 116, 119, 125, 128, 134, 137]
dp_pp.atoms.replace(H_list, H)

dp_pp = dp_pp.tile(2, 1).tile(2, 0)

# dp_pp.write('STRUCT_pp_alter.fdf')

coords = dp_pp.axyz()
ats = list(dp_pp.iter_species())

for ia, a, idx_specie in dp_pp.iter_species():
    if idx_specie == 1:  # H to move
        close = dp_pp.close(ia, R=1.45)
        for ja in close: # C to move towards
            if ats[ja][2] == 0:
                dp_pp = dp_pp.move((1.44 - 1.105) * (coords[ja] - coords[ia]) / 1.44, atoms=ia)

dp_pp = dp_pp.move(-coords[29] + np.array([0.7, 1.0, 10.]))
dp_pp.set_supercell([sc_x, sc_y, sc_z])
dp_pp = remove_outliers(dp_pp)
dp_pp = dp_pp.move(np.array([0., 0., -10.]))
# dp_pp = dp_pp.tile(2, 1).tile(2, 0)

scy_vect = dp_pp.cell[1]
dp_pp = dp_pp.move(-scy_vect, atoms=69).move(-scy_vect, atoms=111)
coords = dp_pp.axyz()

phenyl1 = [13, 15, 14, 16,
           11, 68, 69, 18]
phenyl2 = [71, 72, 21, 22,
           70, 73, 20, 23]
phenyl3 = [106, 108, 45, 47,
           107, 109, 44, 46]
phenyl4 = [53, 55, 50, 52,
           110, 114, 48, 111]

def rotate_phenyl(geom, which, angle):
    '''
    Rotate a phenyl from the planar position, the rotation vectors go from right to left.
    geom must be based on dp_pp_ortho (same indices)
    '''
    assert which in ['1', '2', '3', '4']
    if which == '1':
        rotated = geom.rotate(angle, coords[12] - coords[17], atoms=phenyl1, origo=coords[17], only='xyz')
    elif which == '2':
        rotated = geom.rotate(angle, coords[19] - coords[74], atoms=phenyl2, origo=coords[74], only='xyz')
    elif which == '3':
        rotated = geom.rotate(angle, coords[104] - coords[49], atoms=phenyl3, origo=coords[49], only='xyz')
    else:
        rotated = geom.rotate(angle, coords[51] - coords[54], atoms=phenyl4, origo=coords[54], only='xyz')
    return rotated


def rotate_dppp(geom, ang1, ang2, ang3, ang4):
    rp = rotate_phenyl(geom, '1', ang1)
    rp = rotate_phenyl(rp, '2', ang2)
    rp = rotate_phenyl(rp, '3', ang3)
    return rotate_phenyl(rp, '4', ang4)


# # General creation of structures
# for ang_a in range(-85, 95, 5):
#     for ang_b in range(-85, 95, 5):
#         rotating_dppp = rotate_dppp(dp_pp, ang_a, ang_a + ang_b, -ang_a, -ang_a - ang_b)
#         rotating_dppp.write('structures/STRUCT_{0:03d}_{1:03d}.fdf'.format(ang_a, ang_b))

ang_a = 0
# for j, ang_b in enumerate(range(-75, 105, 15)):
#     rotating_dppp = rotate_dppp(dp_pp, ang_a, ang_a + ang_b, -ang_a, -ang_a - ang_b)
#     rotating_dppp.write('structures/STRUCT_{0:02d}_{1:02d}.fdf'.format(ang_a, j))

rotated = rotate_dppp(dp_pp, 0., 0., 0., 0.)
sisl.plot(rotated, s=10); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# ax.set_xlim(25., 40.)#; ax.set_ylim(0., 10.)
plt.tight_layout(); plt.show()
# sisl.plot(dp_pp, s=5, atom_indices=True); ax = plt.gca(); ax.set_aspect('equal', adjustable='box'); plt.tight_layout()
# fig=plt.gcf(); w, h = fig.get_size_inches(); fig.set_size_inches(w * 2, h * 2); plt.savefig('alter_h.pdf')
