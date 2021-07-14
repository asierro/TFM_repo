import sisl
import numpy as np
import matplotlib.pyplot as plt

'''
Creation of dp_pp structure with different twist angles of the phenyl groups (based off of the relaxed structure)
'''

XV_sile = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/dp_pp/relax/NPG_PP-nonplanar.XV')
dp_pp_ortho = sisl.Geometry.read(XV_sile)

coords = dp_pp_ortho.axyz()

# Move H atoms belonging to phenyls that were on the other side of the supercell
sc_y = dp_pp_ortho.cell[1]
dp_pp_ortho = dp_pp_ortho.move(-sc_y, atoms=19).move(-sc_y, atoms=24)


phenyl1L = [103, 104, 105, 106,
            16, 18, 19, 13]
phenyl2L = [109, 110, 111, 112,
            0, 17, 14, 10]
phenyl1R = [119, 118, 122, 116,
            22, 23, 26, 21]
phenyl2R = [114, 121, 91, 125,
            25, 34, 27, 24]

# rotated1L = dp_pp_ortho.rotate(61.33, coords[102] - coords[107], atoms=phenyl1L, origo=coords[107], only='xyz')
# rotated2L = dp_pp_ortho.rotate(118.23, coords[108] - coords[113], atoms=phenyl2L, origo=coords[113], only='xyz')
# rotated1R = dp_pp_ortho.rotate(61.33, coords[115] - coords[117], atoms=phenyl1R, origo=coords[117], only='xyz')
# rotated2R = dp_pp_ortho.rotate(118.23, coords[90] - coords[120], atoms=phenyl2R, origo=coords[120], only='xyz')

# # sisl.plot(dp_pp_ortho, s=10); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# sisl.plot(rotated1R, s=10); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# plt.tight_layout(); plt.show()


def rotate_phenyl(geom, which, angle):
    '''
    Rotate a phenyl from the planar position, the rotation vectors go from right to left.
    geom must be based on dp_pp_ortho (same indices)
    '''
    assert which in ['1L', '2L', '1R', '2R']
    if which == '1L':
        rotated = geom.rotate(angle - 28.67, coords[102] - coords[107], atoms=phenyl1L, origo=coords[107], only='xyz')
    elif which == '2L':
        rotated = geom.rotate(angle + 28.23, coords[108] - coords[113], atoms=phenyl2L, origo=coords[113], only='xyz')
    elif which == '1R':
        rotated = geom.rotate(angle - 28.28, coords[115] - coords[117], atoms=phenyl1R, origo=coords[117], only='xyz')
    else:
        rotated = geom.rotate(angle + 28.85, coords[90] - coords[120], atoms=phenyl2R, origo=coords[120], only='xyz')
    return rotated


def rotate_dppp(geom, ang1L, ang2L, ang1R, ang2R):
    rp = rotate_phenyl(geom, '1L', ang1L)
    rp = rotate_phenyl(rp, '2L', ang2L)
    rp = rotate_phenyl(rp, '1R', ang1R)
    return rotate_phenyl(rp, '2R', ang2R)

# Create STRUCT files with the rotated dp_pp.
# I will rotate one phenyl aº from the planar position from 0º to 90º,
# and the other bº from aº to 90º+aº.
# The other bridge will be rotated "antisymmetrically" in order for both ribbons to be equivalent.
# Specifically one ribbon is rotated onto the other with a 180º x-axis rotation

for ang_a in range(0, 95, 5):
    for ang_b in range(0, 95, 5):
        rotating_dppp = rotate_dppp(dp_pp_ortho, ang_a, ang_a + ang_b, ang_a + 180, ang_a + ang_b + 180)
        rotating_dppp.write('structures/STRUCT_{}_{}.fdf'.format(ang_a, ang_b))
        # f = plt.figure()
        # f.clear()
        # plt.close(f)
        # plt.clf()
        # sisl.plot(rotating_dppp, s=10); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
        # plt.tight_layout()
        # plt.savefig('/home/asier/Documents/master/tfm/sisl/rotating_dppp/gif/rot_{}_{}.png'.format(ang_a, ang_b))