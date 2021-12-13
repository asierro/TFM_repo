import sisl
import sisl.viz
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from rotating_dnbp2 import torsion_angle

'''
Creation of dp_pp structure functionalized with two nitro groups, replacing the biphenyl bridges with
2,2'-dinitrobiphenyl, with different twist angles of the phenyl groups (while keeping nitro groups planar wrt the 
corresponding ring.
Based off of the relaxed planar structure.
'''


def plotter():
    sisl.plot(dn_dppp, s=10, atom_indices=True)
    ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
    # ax.set_xlim(right=12); ax.set_ylim(top=10)
    plt.show()


def plotlier():
    plot = dn_dppp.plot()
    # plot.update_settings(backend='plotly')
    plot.show()


def torsion_root_left(beta, geom):
    """
    Returns the torsion angle b/n ats 105, 103, 102 and 104 (phenyl 2), depending on the rotation of one nitrophenyl
    (rot_atoms_left) by an angle beta. (This is to find the initial configuration)
    :param geom: geometry of the DN_DPPP (only this specific one works)
    :param beta: angle of rotation of phenyl (in degrees (º))
    :return: rotated geometry
    """
    geom_coords = geom.axyz()
    rot_geom = geom.rotate(beta, geom_coords[102] - geom_coords[117], atoms=rot_atoms_left, origin=geom_coords[117],
                           only='xyz')
    rot_coords = rot_geom.axyz()
    return torsion_angle(rot_coords[105], rot_coords[103], rot_coords[102], rot_coords[104])


def torsion_root_right(beta, geom):
    """
    Returns the torsion angle b/n ats 129, 127, 126 and 128 (phenyl 4), depending on the rotation of one nitrophenyl
    (rot_atoms_left) by an angle beta. (This is to find the initial configuration)
    :param geom: geometry of the DN_DPPP (only this specific one works)
    :param beta: angle of rotation of phenyl (in degrees (º))
    :return: rotated geometry
    """
    geom_coords = geom.axyz()
    rot_geom = geom.rotate(beta, geom_coords[126] - geom_coords[141], atoms=rot_atoms_right, origin=geom_coords[141],
                           only='xyz')
    rot_coords = rot_geom.axyz()
    return torsion_angle(rot_coords[129], rot_coords[127], rot_coords[126], rot_coords[128])


def half_torsion_left1(alpha, geom):
    """
    Returns the torsion angle b/n ats 5, 3, 2 and a point above 2 in the z direction, depending on rotation angle alpha
    about the axis defined by coords[17]-coords[16].
    :param geom: geometry of the DNBP (only this specific one works)
    :param alpha: angle of rotation (in degrees (º))
    :return: rotated geometry
    """
    geom_coords = geom.axyz()
    rot_geom = geom.rotate(alpha, geom_coords[77] - geom_coords[53], atoms=left_bridge, origin=geom_coords[53],
                           only='xyz')
    rot_coords = rot_geom.axyz()
    return torsion_angle(rot_coords[105], rot_coords[103], rot_coords[102], rot_coords[102] + [0., 0., 1.])


def half_torsion_left2(alpha, geom):
    """
    Returns the torsion angle b/n ats 5, 3, 2 and a point above 2 in the z direction, depending on rotation angle alpha
    about the axis defined by coords[17]-coords[16].
    :param geom: geometry of the DNBP (only this specific one works)
    :param alpha: angle of rotation (in degrees (º))
    :return: rotated geometry
    """
    geom_coords = geom.axyz()
    rot_geom = geom.rotate(alpha, geom_coords[77] - geom_coords[53], atoms=left_bridge, origin=geom_coords[53], only='xyz')
    rot_coords = rot_geom.axyz()
    return torsion_angle(rot_coords[103] + [0., 0., 1.], rot_coords[103], rot_coords[102], rot_coords[104])


def torsion_diff_left(alpha, geom):
    return half_torsion_left2(alpha, geom) - half_torsion_left1(alpha, geom)


def symmetrize_left(geom):
    """
    Rotate the DNBP geometry wrt axis defined by outer carbons, st the nitro groups lie symmetrically wrt E field (z
    axis). For that we try to find the angle of rotation at which halftorsion1=halftorsion2.
    :param geom: DNBP geometry (specifically this one)
    :return: symmetrized DNBP geometry
    """
    geom_coords = geom.axyz()

    root_results = optimize.root_scalar(torsion_diff_left, args=(geom,), x0=0., x1=10., xtol=0.1)
    if not root_results.converged:
        raise Exception("Root finder for symmetrization did not converge")
    sol = root_results.root

    return geom.rotate(sol, geom_coords[77] - geom_coords[53], atoms=left_bridge, origin=geom_coords[53], only='xyz')


def half_torsion_right1(alpha, geom):
    """
    Returns the torsion angle b/n ats 5, 3, 2 and a point above 2 in the z direction, depending on rotation angle alpha
    about the axis defined by coords[17]-coords[16].
    :param geom: geometry of the DNBP (only this specific one works)
    :param alpha: angle of rotation (in degrees (º))
    :return: rotated geometry
    """
    geom_coords = geom.axyz()
    rot_geom = geom.rotate(alpha, geom_coords[55] - geom_coords[85], atoms=right_bridge, origin=geom_coords[85],
                           only='xyz')
    rot_coords = rot_geom.axyz()
    return torsion_angle(rot_coords[129], rot_coords[127], rot_coords[126], rot_coords[126] + [0., 0., 1.])


def half_torsion_right2(alpha, geom):
    """
    Returns the torsion angle b/n ats 5, 3, 2 and a point above 2 in the z direction, depending on rotation angle alpha
    about the axis defined by coords[17]-coords[16].
    :param geom: geometry of the DNBP (only this specific one works)
    :param alpha: angle of rotation (in degrees (º))
    :return: rotated geometry
    """
    geom_coords = geom.axyz()
    rot_geom = geom.rotate(alpha, geom_coords[55] - geom_coords[85], atoms=right_bridge, origin=geom_coords[85], only='xyz')
    rot_coords = rot_geom.axyz()
    return torsion_angle(rot_coords[127] + [0., 0., 1.], rot_coords[127], rot_coords[126], rot_coords[128])


def torsion_diff_right(alpha, geom):
    return half_torsion_right2(alpha, geom) - half_torsion_right1(alpha, geom)


def symmetrize_right(geom):
    """
    Rotate the DNBP geometry wrt axis defined by outer carbons, st the nitro groups lie symmetrically wrt E field (z
    axis). For that we try to find the angle of rotation at which halftorsion1=halftorsion2.
    :param geom: DNBP geometry (specifically this one)
    :return: symmetrized DNBP geometry
    """
    geom_coords = geom.axyz()

    root_results = optimize.root_scalar(torsion_diff_right, args=(geom,), x0=0., x1=10., xtol=0.1)
    if not root_results.converged:
        raise Exception("Root finder for symmetrization did not converge")
    sol = root_results.root

    # Add 180 so it is upside down wrt the left bridge
    return geom.rotate(sol+180, geom_coords[55] - geom_coords[85], atoms=right_bridge, origin=geom_coords[85], only='xyz')


def too_close(coordinates):
    for i in nitro1:
        for j in nitro2:
            if np.linalg.norm(coordinates[i] - coordinates[j]) < 1.4:  # Nitro groups too close
                return True
    for i in nitro3:
        for j in nitro4:
            if np.linalg.norm(coordinates[i] - coordinates[j]) < 1.4:  # Nitro groups too close
                return True
    if np.linalg.norm(coordinates[109] - coordinates[118]) < 1.1 or \
       np.linalg.norm(coordinates[108] - coordinates[119]) < 1.1 or \
       np.linalg.norm(coordinates[133] - coordinates[142]) < 1.1 or \
       np.linalg.norm(coordinates[132] - coordinates[143]) < 1.1:  # oxygen too close to hydrogen atom
        return True
    return False


XV_sile = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/dn_dppp/DN_DPPP.XV')
dn_dppp = sisl.Geometry.read(XV_sile)
coords = dn_dppp.axyz()

# torsion_atoms_left = [105, 103, 102, 104]
# torsion_atoms_right = [129, 127, 126, 128]
left_bridge = [116, 113, 114, 105, 103, 111,
               119, 122, 121, 100, 107, 109,
               102, 104, 110, 112, 115, 117,
               118, 123, 120, 101, 108, 106]
right_bridge = [127, 129, 135, 137, 138, 140,
                145, 146, 143, 124, 131, 133,
                126, 128, 134, 136, 139, 141,
                144, 147, 142, 125, 132, 130]

rot_atoms_left = [110, 115, 104, 112,
                  118, 123, 120,
                  101, 108, 106]
rot_atoms_right = [128, 136, 139, 134,
                   142, 147, 144,
                   125, 132, 130]
nitro1 = [100, 107, 109]; nitro2 = [101, 108, 106]
nitro3 = [124, 131, 133]; nitro4 = [125, 132, 130]

# Rotate phenyl 2 st the initial configuration in the left bridge
# has the two nitros lined up (this way, a rotation of 180º will be
# enough to cover all possible configurations
left_root_results = optimize.root_scalar(torsion_root_left, args=(dn_dppp,), x0=-90., x1=-100., xtol=0.1)
beta_left = left_root_results.root
dn_dppp = dn_dppp.rotate(beta_left, coords[102] - coords[117], atoms=rot_atoms_left, origin=coords[117], only='xyz')

coords = dn_dppp.axyz()

right_root_results = optimize.root_scalar(torsion_root_right, args=(dn_dppp,), x0=-90., x1=-100., xtol=0.1)
beta_right = right_root_results.root
dn_dppp = dn_dppp.rotate(beta_right, coords[126] - coords[141], atoms=rot_atoms_right, origin=coords[141], only='xyz')

dn_dppp = symmetrize_left(dn_dppp)
dn_dppp = symmetrize_right(dn_dppp)

coords = dn_dppp.axyz()

for ang in range(0, 190, 10):
    rotated = dn_dppp.rotate(ang, coords[102] - coords[117], atoms=rot_atoms_left, origin=coords[117], only='xyz')
    rotated = symmetrize_left(rotated)
    rotated = symmetrize_left(rotated)
    rotated = rotated.rotate(ang, coords[126] - coords[141], atoms=rot_atoms_right, origin=coords[141], only='xyz')
    rotated = symmetrize_right(rotated)
    rotated = symmetrize_right(rotated)

    new_coords = rotated.axyz()
    if too_close(new_coords):
        # print('Too close!')
        continue
    print(str(ang) + ' ', torsion_angle(new_coords[105], new_coords[103], new_coords[102], new_coords[104]), ' ',
          torsion_angle(new_coords[129], new_coords[127], new_coords[126], new_coords[128]))

    # rotated.write('structures/STRUCT_{0:03d}.fdf'.format(ang))

