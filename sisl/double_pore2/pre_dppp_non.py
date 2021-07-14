import sisl
import numpy as np
import matplotlib.pyplot as plt
'''
Prepare a non orthogonal DP PP for relaxation, based on the relaxed orthogonal DP PP
'''

XV_sile = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/dp_pp/relax/NPG_PP-nonplanar.XV')

dp_pp_ortho = sisl.Geometry.read(XV_sile)

dp_pp_ortho = dp_pp_ortho.tile(2, 1)

c = dp_pp_ortho.axyz()

dp_pp_ortho = dp_pp_ortho.move(-c[63] + np.array([0., 0., 25.])) #, cell=True) #move up 25A for atoms to be inside sc

c = dp_pp_ortho.axyz()

sc_1 = c[203]; sc_1[2]=0.
sc_2 = c[80]; sc_2[2]=0.
sc_3 = [0., 0., 50.]
dp_pp_ortho.set_supercell([sc_1, sc_2, sc_3])

c = dp_pp_ortho.axyz()


# def outside_sc(point, sc_1, sc_2, sc_3):
#     u, v, w = (np.cross(sc_1, sc_2), np.cross(sc_1, sc_3), np.cross(sc_2, sc_3))
#     if (0 <= np.dot(u, point) < np.dot(u, sc_3) and
#             0 <= np.dot(v, point) < np.dot(v, sc_2) and
#             0 <= np.dot(w, point) < np.dot(w, sc_1)):
#         return True
#     return False


def remove_outliers(geom):
    '''
    This routine removes all atoms in a geometry which are located outside the supercell (the parallepiped defined by
    the sc vectors)
    '''
    sc = geom.cell
    print('sc=',sc)
    axyz = geom.axyz()
    # deal with the signs in a more systematic way
    u, v, w = (-np.cross(sc[0], sc[1]), np.cross(sc[0], sc[2]), -np.cross(sc[1], sc[2]))
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


dp_pp_non = remove_outliers(dp_pp_ortho)
#dp_pp_non = dp_pp_non.tile(2, 1).tile(2, 0)

sisl.plot(dp_pp_non, s=10); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
plt.tight_layout(); plt.show()
