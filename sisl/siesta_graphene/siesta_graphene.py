import numpy as np
from sisl import *
import matplotlib.pyplot as plt

graphene = geom.graphene(1.44)
# plot(graphene)
# plt.show()

# open('RUN.fdf', 'w').write("""%include STRUCT.fdf
# SystemLabel siesta_2
# PAO.BasisSize SZP
# MeshCutoff 250. Ry
# CDF.Save true
# CDF.Compress 9
# SaveHS true
# SaveRho true
# %block kgrid.MonkhorstPack
#   61  1 1 0.
#    1 61 1 0.
#    0  0 1 0.
# %endblock
# """)
# graphene.write('STRUCT.fdf')

fdf = get_sile('RUN.fdf')
H = fdf.read_hamiltonian()
# print(H)

# E = np.linspace(-6, 4, 500)
# # for nk in [21, 41, 61, 81]:
# #     bz = MonkhorstPack(H, [nk, nk, 1])
# #     plt.plot(E, bz.apply.average.DOS(E), label='nk={}'.format(nk));
# # plt.xlim(E[0], E[-1]); plt.ylim(0, None)
# # plt.xlabel(r'$E - E_F$ [eV]')
# # plt.ylabel(r'DOS [1/eV]')
# # plt.legend();
# # plt.show()
#
# bz = MonkhorstPack(H, [81, 81, 1])
# idx_s = list()
# idx_pxy = list()
# idx_pz = list()
# for i, orb in enumerate(H.geometry.atoms[0]):
#     if orb.l == 0:
#         idx_s.append(i)
#     elif orb.l == 1 and (orb.m in [-1, 1]):
#         idx_pxy.append(i)
#     elif orb.l == 1 and orb.m == 0:
#         idx_pz.append(i)
# print('Orbital index of s: {}'.format(idx_s))
# print('Orbital index of p_x and p_y: {}'.format(idx_pxy))
# print('Orbital index of p_z: {}'.format(idx_pz))
# # Get all orbitals
# all_s = np.add.outer(H.geometry.a2o([0, 1]), idx_s).ravel()
# all_pxy = np.add.outer(H.geometry.a2o([0, 1]), idx_pxy).ravel()
# all_pz = np.add.outer(H.geometry.a2o([0, 1]), idx_pz).ravel()
# def wrap(PDOS):
#     pdos_s = PDOS[all_s, :].sum(0)
#     pdos_pxy = PDOS[all_pxy, :].sum(0)
#     pdos_pz = PDOS[all_pz, :].sum(0)
#     return np.stack((pdos_s, pdos_pxy, pdos_pz))
# # pDOS = bz.apply.average.PDOS(E, wrap=wrap)
# # plt.plot(E, pDOS[0, :], label='$s$');
# # plt.plot(E, pDOS[1, :], label='$p_x+p_y$');
# # plt.plot(E, pDOS[2, :], label=r'$p_z$');
# # plt.xlim(E[0], E[-1]); plt.ylim(0, None)
# # plt.xlabel(r'$E - E_F$ [eV]')
# # plt.ylabel(r'DOS [1/eV]')
# # plt.legend();
# # plt.show()

weight_s = list()
weight_pxy = list()
weight_pz = list()
def wrap_fatbands(eigenstate):
    # The eigenstate object has several features.
    # For now we will simply calculate the <psi_i| S(k) | psi_i> weight for
    # the orbitals we are interested in.
    norm2 = eigenstate.norm2(sum=False)
    weight_s.append(norm2[:, all_s].sum(-1))
    weight_pxy.append(norm2[:, all_pxy].sum(-1))
    weight_pz.append(norm2[:, all_pz].sum(-1))
    return eigenstate.eig
# Define the band-structure
bz = BandStructure(H, [[0] * 3, [2./3, 1./3, 0], [0.5, 0.5, 0], [1] * 3], 400,
                   name=[r'$\Gamma$', r'$K$', r'$M$', r'$\Gamma$'])

onlybands = True

if not onlybands:
    # Calculate all eigenvalues
    eig = bz.apply.array.eigenstate(wrap=wrap_fatbands).T
    # eig = bz.apply.array.eigh().T

    weight_s = np.array(weight_s).T
    weight_pxy = np.array(weight_pxy).T
    weight_pz = np.array(weight_pz).T

    linear_k, k_tick, k_label = bz.lineark(True)

    Emin, Emax = -21, 10
    # This is to determine the width of the fat-bands
    # The width of the fat-bands is dependent on the energy range and also on the variety
    # of contributions.
    dE = (Emax - Emin) / 20.
    plt.ylabel(r'$E-E_F$ [eV]')
    plt.xlim(linear_k[0], linear_k[-1])
    plt.xticks(k_tick, k_label)
    plt.ylim(Emin, Emax)

    # Now plot the bands
    for i, e in enumerate(eig):
        s = np.abs(weight_s[i, :] * dE)
        pxy = np.abs(weight_pxy[i, :] * dE)
        pz = np.abs(weight_pz[i, :] * dE)
        plt.plot(linear_k, e, color='k') # black-line (band-structure)
        # Full fat-band
        plt.fill_between(linear_k, e - dE, e + dE, color='k', alpha=0.1)
        # pz
        plt.fill_between(linear_k, e + (s + pxy), e + (s + pxy + pz), color='g', alpha=0.5)
        plt.fill_between(linear_k, e - (s + pxy + pz), e - (s + pxy), color='g', alpha=0.5)
        # pxy
        plt.fill_between(linear_k, e + (s), e + (s + pxy), color='b', alpha=0.5)
        plt.fill_between(linear_k, e - (s), e - (s + pxy), color='b', alpha=0.5)
        # s
        plt.fill_between(linear_k, e - (s), e + (s), color='r', alpha=0.5)
    plt.show()
else:
    eig = bz.apply.array.eigh().T

    linear_k, k_tick, k_label = bz.lineark(True)

    Emin, Emax = -21, 10
    plt.ylabel(r'$E-E_F$ [eV]')
    plt.xlim(linear_k[0], linear_k[-1])
    plt.xticks(k_tick, k_label)
    plt.ylim(Emin, Emax)

    # Now plot the bands
    for e in eig:
        plt.plot(linear_k, e, color='k')
    plt.show()
