import sisl
import numpy as np
import matplotlib.pyplot as plt

graphene = sisl.geom.graphene()
# print(graphene)
# sisl.plot(graphene)
# plt.show()

H = sisl.Hamiltonian(graphene)
# print(H)

# # First method to specify matrix elements
# H[0, 0] = 0.
# H[1, 1] = 0.
# H[0, 1] = -2.7
# H[1, 0] = -2.7
#
# H[0, 1, (-1, 0)] = -2.7
# H[0, 1, (0, -1)] = -2.7
# H[1, 0, (1, 0)] = -2.7
# H[1, 0, (0, 1)] = -2.7

# # 2nd method
# for ia, io in H.iter_orbitals():
#     print('a')
#     idx = H.geom.close(ia, R=[0.1, 1.43])
#     print(io, idx[0])
#     print(io, idx[1])
#     H[io, idx[0]] = 0.
#     H[io, idx[1]] = -2.7

# 3rd method
H.construct([[0.1, 1.43], [0., -2.7]])

# print("Gamma point:", H.eigh())
# print("K point:", H.eigh(k=[2./3,1./3,0]))

band = sisl.BandStructure(H, [[0., 0.], [2./3, 1./3],
                              [1./2, 1./2], [1., 1.]], 301,
                              [r'$\Gamma$', 'K', 'M', r'$\Gamma$'])
eigs = band.apply.array.eigh()

# Retrieve the tick-marks and the linear k points
xtick, xtick_label = band.lineartick()
lk = band.lineark()
plt.plot(lk, eigs)

plt.ylabel('Eigenspectrum [eV]')
plt.gca().xaxis.set_ticks(xtick)
plt.gca().set_xticklabels(xtick_label)

# Also plot x-major lines at the ticks
ymin, ymax = plt.gca().get_ylim()
for tick in xtick:
    plt.plot([tick,tick], [ymin,ymax], 'k')

plt.show()
