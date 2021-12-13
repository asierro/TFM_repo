import sisl
import sisl.viz
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from itertools import chain

# THIS IS HERE JUST FOR IMPORTING FUNCTIONS

'''
Rotate DNBP to different torsion angles and create SIESTA input for relaxation with the constraint of the torsional
angle b/n phenyls (via the ZMatrix), for all unique angles except for those which leave the nitro groups too close 
together.
THIS VERSION IS JUST TO FIND OUT THE ROTATION ANGLES TO SYMMETRIZE DNBP FOR THE DIPOLE APPROXIMATION
'''


def plotter():
    sisl.plot(dnbp, atom_indices=True); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(right=12); ax.set_ylim(top=10)
    plt.show()


def plotlier():
    plot = dnbp.plot()
    plot.show()


def torsion_angle(p0, p1, p2, p3):
    """
    Calculates torsion angle. Source:
    stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    """
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b1 /= np.linalg.norm(b1)

    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def too_close(coordinates):
    for i in nitro1:
        for j in nitro2:
            if np.linalg.norm(coordinates[i] - coordinates[j]) < 1.4:  # Nitro groups too close
                return True
    if np.linalg.norm(coordinates[9] - coordinates[18]) < 1.1 or \
       np.linalg.norm(coordinates[8] - coordinates[19]) < 1.1:  # nitrogen too close to hydrogen atom
        return True
    return False


def half_torsion1(alpha, geom):
    """
    Returns the torsion angle b/n ats 5, 3, 2 and a point above 2 in the z direction, depending on rotation angle alpha
    about the axis defined by coords[17]-coords[16].
    :param geom: geometry of the DNBP (only this specific one works)
    :param alpha: angle of rotation (in degrees (º))
    :return: rotated geometry
    """
    geom_coords = geom.axyz()
    rot_geom = geom.rotate(alpha, geom_coords[17] - geom_coords[16], origo=geom_coords[16], only='xyz')
    rot_coords = rot_geom.axyz()
    return torsion_angle(rot_coords[5], rot_coords[3], rot_coords[2], rot_coords[2] + [0., 0., 1.])


def half_torsion2(alpha, geom):
    """
    Returns the torsion angle b/n ats 5, 3, 2 and a point above 2 in the z direction, depending on rotation angle alpha
    about the axis defined by coords[17]-coords[16].
    :param geom: geometry of the DNBP (only this specific one works)
    :param alpha: angle of rotation (in degrees (º))
    :return: rotated geometry
    """
    geom_coords = geom.axyz()
    rot_geom = geom.rotate(alpha, geom_coords[17] - geom_coords[16], origo=geom_coords[16], only='xyz')
    rot_coords = rot_geom.axyz()
    return torsion_angle(rot_coords[3] + [0., 0., 1.], rot_coords[3], rot_coords[2], rot_coords[4])


# One can check that half_torsion1+halftorsion2(-360º) = torsion angle


def torsion_diff(alpha, geom):
    return half_torsion2(alpha, geom) - half_torsion1(alpha, geom)


def symmetrize(geom):
    """
    Rotate the DNBP geometry wrt axis defined by outer carbons, st the nitro groups lie symmetrically wrt E field (z
    axis). For that we try to find the angle of rotation at which halftorsion1=halftorsion2.
    THIS TIME FORCE THE MOLECULE TO HAVE THE NITROS POINTING UP
    :param geom: DNBP geometry (specifically this one)
    :return: symmetrized DNBP geometry
    """
    geom_coords = geom.axyz()

    root_results = optimize.root_scalar(torsion_diff, args=(geom,), x0=0., x1=10., xtol=0.1)
    if not root_results.converged:
        raise Exception("Root finder for symmetrization did not converge")
    sol = root_results.root
    rotated = geom.rotate(sol, geom_coords[17] - geom_coords[16], origo=geom_coords[16], only='xyz')
    rot_coords = rotated.axyz()
    if rot_coords[5, 2] < rot_coords[3, 2]:
        return rotated.rotate(180, geom_coords[17] - geom_coords[16], origo=geom_coords[16], only='xyz'), sol+180
    return rotated, sol


def torsion_root(beta, geom):
    """
    Returns the torsion angle b/n ats 5, 3, 2 and 4, depending on the rotation of one nitrophenyl (rot_atoms) by an
    angle beta. (This is to find the initial configuration)
    :param geom: geometry of the DNBP (only this specific one works)
    :param beta: angle of rotation of phenyl (in degrees (º))
    :return: rotated geometry
    """
    geom_coords = geom.axyz()
    rot_geom = geom.rotate(beta, geom_coords[2] - geom_coords[17], atoms=rot_atoms, origo=geom_coords[17], only='xyz')
    rot_coords = rot_geom.axyz()
    return torsion_angle(rot_coords[5], rot_coords[3], rot_coords[2], rot_coords[4])


# dbdt_65 is actually DNBP (without thiols)
fdf = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/dbdt_65/RUN.fdf')
dnbp = sisl.Geometry.read(fdf, output=True)  # read from output files

rot_atoms = [4, 12, 10, 15, 18, 23, 20, 1, 8, 6]
all_atoms = list(range(26))
nitro1 = [7, 0, 9]
nitro2 = [8, 1, 6]

coords = dnbp.axyz()

# Rotate the whole DNBP st outermost carbons have the same z coord (so we can apply E-field correctly)
vec = coords[17] - coords[16]
ang = np.arctan2(vec[2], np.linalg.norm(vec[:2]))
rot_vec = [-vec[1], vec[0], 0.]
dnbp = dnbp.rotate(ang, rot_vec, origo=coords[16], only='xyz', rad=True)

dnbp, foo = symmetrize(dnbp)
dnbp, foo = symmetrize(dnbp)
coords = dnbp.axyz()
print('_________________________')

# WE ARE NOT DOING THIS NOW BC WE DIDN'T DO IT IN THE OLD EFIELD=0 CASE
# # Rotate nitrophenyl st the initial configuration has the two nitros lined up (this way, a rotation of 180º will be
# # enough to cover all possible configurations
# beta_root_results = optimize.root_scalar(torsion_root, args=(dnbp,), x0=-110., x1=-115., xtol=0.1)
# beta = beta_root_results.root
# dnbp = dnbp.rotate(beta, coords[2] - coords[17], atoms=rot_atoms, origin=coords[17], only='xyz')
#
# dnbp = symmetrize(dnbp)
# dnbp = symmetrize(dnbp)
# coords = dnbp.axyz()

for angle in range(0, 360, 5):
    rotated = dnbp.rotate(angle, coords[2] - coords[17], atoms=rot_atoms, origo=coords[17], only='xyz')
    new_coords = rotated.axyz()
    if too_close(new_coords):
        # print('Too close!')
        continue
    rotated, sol1 = symmetrize(rotated)
    rotated, sol2 = symmetrize(rotated)
    new_coords = rotated.axyz()
    print(str(angle) + ' ', torsion_angle(new_coords[5], new_coords[3], new_coords[2], new_coords[4]), ' ', sol1 + sol2)

# sol1+sol2 is the angle which the molecule needs to be rotated by for it to be symmetrized. This is for the dipole
# approximation (to rotate the dipoles in the same way)

