import sisl
import sisl.viz
import numpy as np
import matplotlib.pyplot as plt

'''
Slater-Koster tight-binding Hamiltonian to try to explain 75º minimum of delta E at gamma.
'''

XV_sile = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/dp_pp/relax/NPG_PP-nonplanar.XV')
dppp = sisl.Geometry.read(XV_sile)

# Move H atoms belonging to phenyls that were on the other side of the supercell
sc_y = dppp.cell[1]
dppp = dppp.move(-sc_y, atoms=19).move(-sc_y, atoms=24)
coords = dppp.axyz()

# Remove all H atoms
H_remove = []
for ia, a, idx_specie in dppp.iter_species():
    if idx_specie == 0:
        H_remove.append(ia)

dppp = dppp.remove(H_remove)
coords = dppp.axyz()

rot_phenyl1 = [67, 68, 69, 70]
rot_phenyl2 = [73, 74, 75, 76]
rot_phenyl3 = [83, 82, 86, 80]
rot_phenyl4 = [78, 85, 55, 89]

# Phenyls including the ones on each end:
phenyl1 = [67, 68, 69, 70, 66, 71]
phenyl2 = [73, 74, 75, 76, 72, 77]
phenyl3 = [83, 82, 86, 80, 79, 81]
phenyl4 = [78, 85, 55, 89, 54, 84]
p24 = phenyl2 + phenyl4

# plot = dppp.plot()
# plot.show()


def rotate_phenyl(geom, which, angle):
    '''
    Rotate a phenyl from the planar position, the rotation vectors go from right to left.
    geom must be based on dp_pp_ortho (same indices)
    '''
    assert which in [1, 2, 3, 4]
    if which == 1:
        rotated = geom.rotate(angle - 28.67, coords[66] - coords[71], atoms=rot_phenyl1, origin=coords[71], only='xyz')
    elif which == 2:
        rotated = geom.rotate(angle + 28.23, coords[72] - coords[77], atoms=rot_phenyl2, origin=coords[77], only='xyz')
    elif which == 3:
        rotated = geom.rotate(angle - 28.28, coords[79] - coords[81], atoms=rot_phenyl3, origin=coords[81], only='xyz')
    else:
        rotated = geom.rotate(angle + 28.85, coords[54] - coords[84], atoms=rot_phenyl4, origin=coords[84], only='xyz')
    return rotated


def rotate_dppp(geom, ang1, ang2, ang3, ang4):
    rp = rotate_phenyl(geom, 1, ang1)
    rp = rotate_phenyl(rp, 2, ang2)
    rp = rotate_phenyl(rp, 3, ang3)
    return rotate_phenyl(rp, 4, ang4)


def Exx(lmn, Vpp_pi, Vpp_sigma):
    l = lmn[0]
    return l**2 * Vpp_sigma + (1 - l**2) * Vpp_pi


def Exy(lmn, Vpp_pi, Vpp_sigma):
    l = lmn[0]
    m = lmn[1]
    return l * m * (Vpp_sigma - Vpp_pi)


def Exz(lmn, Vpp_pi, Vpp_sigma):
    l = lmn[0]
    n = lmn[2]
    return l * n * (Vpp_sigma - Vpp_pi)


def Eij(lmn, pz_spherical, Vpp_pi, Vpp_sigma):
    theta = pz_spherical[1]
    phi = pz_spherical[2]
    return (Exy(lmn, Vpp_pi, Vpp_sigma) * np.cos(phi) * np.sin(theta) +
            Exz(lmn, Vpp_pi, Vpp_sigma) * np.sin(phi) * np.sin(theta) +
            Exx(lmn, Vpp_pi, Vpp_sigma) * np.cos(theta))


def to_spherical(v):
    xy = v[0]**2 + v[1]**2
    return np.array([np.sqrt(xy + v[2] ** 2),  # r
                     np.arctan2(np.sqrt(xy), v[2]),  # theta
                     np.arctan2(v[1], v[0])])  # phi


def rodrigues(v, k, alpha):
    '''
    Return the rotated vector v about vector k by an angle alpha, using the Rodrigues formula
    '''
    k = k / np.linalg.norm(k)
    return v * np.cos(alpha) + np.cross(k, v) * np.sin(alpha) + k * k.dot(v) * (1 - np.cos(alpha))

# ---------------------------------------------------------------------------------------------------------------------
# Create TB Hamiltonian (rotate beforehand, this is just for testing)
# MaxR is negative for some reason so replace the C Atoms in this geometry for C Atoms with the correct orbital radius
# I'm going to set an e^-r radial dependence of the orbital, even though this does not actually matter because we are
# going to parametrize hoppings
# maxR = 3.1
maxR = 4.6
# r = np.linspace(0, maxR, 50)
# f = np.exp(-r)
# orb = sisl.AtomicOrbital('2pz', (r, f))
# C = sisl.Atom(Z=6, orbitals=orb)
C = sisl.Atom(Z=6, R=maxR)
atoms = sisl.Atoms(atoms=[C], na=104)
dppp = sisl.Geometry(coords, atoms=atoms, sc=dppp.cell)
dppp.set_nsc([3, 3, 1])

# "Flat" DPPP for when both atoms being considered are in the rotating phenyl
flat_dppp = rotate_dppp(dppp, 0., 0., 180., 180.)
# flat_dppp.set_nsc([3, 3, 1])

ang_a = 0.
for ang_b in range(0, 95, 5):
    rotating_dppp = rotate_dppp(dppp, ang_a, ang_a + ang_b, ang_a + 180, ang_a + ang_b + 180)
    plot = rotating_dppp.plot()
    plot.show()

    H_SK = sisl.Hamiltonian(rotating_dppp)

    z = np.array([0, 0, 1])

    # Direction vectors of pz orbitals in phenyls 2 and 4 after rotating them
    pz_vect_p2 = to_spherical(rodrigues(z, coords[72] - coords[77], ang_b*np.pi/180))
    pz_vect_p4 = to_spherical(rodrigues(z, coords[54] - coords[84], (ang_b + 180.)*np.pi/180))

    a = 1.45  # avg lattice constant
    d_2nn = 2 * np.cos(30 * np.pi / 180) * a  # ~dist b/n 2nd neighbours
    t_pi = -2.7  # 2.682
    t_sigma = -.371 # .8726  # find better params?
    t_2nn = -.0027  # 0.2  # 2NN hopping
    q_pi = np.log(t_2nn / t_pi) / (1 - d_2nn / a)  # decay rate. Going to use q_sigma = q_pi for now
    q_sigma = 3.34 / 1.42 * q_pi


    for ia in rotating_dppp:
        idx, xyz, rij = rotating_dppp.close(ia, R=(0.1, maxR + 0.2), ret_xyz=True, ret_rij=True)
        # idx_flat, xyz_flat = flat_dppp.close(ia, R=(0.1, maxR + 0.2), ret_xyz=True)
        # idx: list with two els: 1: onsite index (ia), 2: indices of neighbours w/in maxR+0.1
        # xyz: list with two els: 1: ia atom pos, 2: pos of neighbours
        # xyz: list with two els: 1: dist from ia to ia (0.), 2: dist to neighbours
        H_SK[ia, idx[0]] = 0.

        vij, rij = xyz[1] - rotating_dppp.xyz[ia], rij[1]  # relative pos of neighbours wrt ia
        # vij_flat = xyz_flat[1] - flat_dppp.xyz[ia]

        # Calculate direction cosines:
        l = vij.dot(np.array([0., 0., 1.])) / np.sqrt((vij ** 2).sum(axis=1))
        m = vij.dot(np.array([1., 0., 0.])) / np.sqrt((vij ** 2).sum(axis=1))
        n = vij.dot(np.array([0., 1., 0.])) / np.sqrt((vij ** 2).sum(axis=1))

        # l_flat = vij_flat.dot(np.array([0., 0., 1.])) / np.sqrt((vij_flat ** 2).sum(axis=1))
        # m_flat = vij_flat.dot(np.array([1., 0., 0.])) / np.sqrt((vij_flat ** 2).sum(axis=1))
        # n_flat = vij_flat.dot(np.array([0., 1., 0.])) / np.sqrt((vij_flat ** 2).sum(axis=1))

        Vpp_pi = t_pi * np.exp(q_pi * (1 - rij / a))
        Vpp_sigma = t_sigma * np.exp(q_sigma * (1 - rij / a))

        for j, ja in enumerate(idx[1]):
            lmn = [l[j], m[j], n[j]]
            # lmn_flat = [l_flat[j], m_flat[j], n_flat[j]]
            Vpi_j = Vpp_pi[j]
            Vsigma_j = Vpp_sigma[j]
            if ia not in p24 and ja in p24:
                if ja in phenyl2:
                    H_SK[ia, ja] = Eij(lmn, pz_vect_p2, Vpi_j, Vsigma_j)
                elif ja in phenyl4:
                    H_SK[ia, ja] = Eij(lmn, pz_vect_p4, Vpi_j, Vsigma_j)
            elif ia in p24 and ja not in p24:
                if ia in phenyl2:
                    pz_vect = pz_vect_p2 + np.array([0., 0., np.pi])
                    H_SK[ia, ja] = Eij(lmn, pz_vect, Vpi_j, Vsigma_j)
                elif ia in phenyl4:
                    pz_vect = pz_vect_p4 + np.array([0., 0., np.pi])
                    H_SK[ia, ja] = Eij(lmn, pz_vect, Vpi_j, Vsigma_j)
            elif ia in p24 and ja in p24:
                # Bc of the max radius this only happens when they are in the same phenyl (2 or 4)
                xyz_ia_flat = flat_dppp.axyz()[ia]
                xyz_ja_flat = flat_dppp.axyz()[ja]
                vij_flat = xyz_ja_flat - xyz_ia_flat
                l_flat = vij_flat.dot(np.array([0., 0., 1.])) / np.sqrt((vij_flat ** 2).sum())
                m_flat = vij_flat.dot(np.array([1., 0., 0.])) / np.sqrt((vij_flat ** 2).sum())
                n_flat = vij_flat.dot(np.array([0., 1., 0.])) / np.sqrt((vij_flat ** 2).sum())
                lmn_flat = [l_flat, m_flat, n_flat]
                if ja in phenyl2:
                    H_SK[ia, ja] = Exx(lmn_flat, Vpi_j, Vsigma_j)
                elif ja in phenyl4:
                    H_SK[ia, ja] = Exx(lmn_flat, Vpi_j, Vsigma_j)
            else:
                H_SK[ia, ja] = Exx(lmn, Vpi_j, Vsigma_j)


    # Band structure plotting
    band = sisl.BandStructure(H_SK, [[0.5, 0., 0.], [0., 0., 0.],[0., 0.5, 0.]], 301, ['X', r'$\Gamma$', 'Y'])
    eigs = band.apply.array.eigh(gauge='r')

    print(eigs[46, 53] - eigs[46, 52])

# # Retrieve the tick-marks and the linear k points
# xtick, xtick_label = band.lineartick()
# lk = band.lineark()
# plt.plot(lk, eigs)
#
# plt.ylabel('Eigenspectrum [eV]')
# plt.gca().set_ylim(0.5, 0.7)#(-1., 3.)  # (0.5, 0.7)
# plt.gca().xaxis.set_ticks(xtick)
# plt.gca().set_xticklabels(xtick_label)
#
# # Also plot x-major lines at the ticks
# ymin, ymax = plt.gca().get_ylim()
# for tick in xtick:
#     plt.plot([tick, tick], [ymin, ymax], 'k')
#
# plt.show()

# I will rotate one phenyl aº from the planar position from 0º to 90º,
# and the other bº from aº to 90º+aº.
# The other bridge will be rotated "antisymmetrically" in order for both ribbons to be equivalent.
# Specifically one ribbon is rotated onto the other with a 180º x-axis rotation

# for ang_a in range(0, 95, 5):
#     for ang_b in range(0, 95, 5):
#         rotating_dppp = rotate_dppp(dp_pp_ortho, ang_a, ang_a + ang_b, ang_a + 180, ang_a + ang_b + 180)
#         rotating_dppp.write('structures/STRUCT_{}_{}.fdf'.format(ang_a, ang_b))
        # f = plt.figure()
        # f.clear()
        # plt.close(f)
        # plt.clf()
        # sisl.plot(rotating_dppp, s=10); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
        # plt.tight_layout()
        # plt.savefig('/home/asier/Documents/master/tfm/sisl/rotating_dppp/gif/rot_{}_{}.png'.format(ang_a, ang_b))