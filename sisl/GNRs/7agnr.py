import sisl
import numpy as np
import matplotlib.pyplot as plt

C = sisl.Atom(6, mass=12.0107, tag='C')
agnr7 = sisl.geom.nanoribbon(bond=1.43, atoms=C, width=7, kind='armchair')

# plt.gca().set_aspect('equal', adjustable='box')
# sisl.plot(agnr7); plt.show()

# open('RUN.fdf', 'w').write("""%include STRUCT.fdf
# SystemLabel agnr7
# PAO.BasisSize SZP
# MeshCutoff 250. Ry
# CDF.Save true
# CDF.Compress 9
# SaveHS true
# SaveRho true
# %block kgrid.MonkhorstPack
#   101  1 1 0.
#    0 1 0 0.
#    0  0 1 0.
# %endblock
# """)
# agnr7.write('STRUCT.fdf')

fdf = sisl.get_sile('RUN.fdf')
H = fdf.read_hamiltonian()

band = sisl.BandStructure(H,
    [[-1/2, 0, 0],
     [0, 0, 0],
     [1/2, 0, 0]],
    400)

# lk, kt, kl = band.lineark(True)
# bs = band.apply.array.eigh(eta=True)
# plt.plot(lk, bs, 'k', label=r'SZP')
# plt.ylabel(r'$E-E_F$ [eV]'); plt.xticks(kt, kl)
# for x in kt:
#     plt.axvline(x, c='k', ls=':', lw=0.5)
# plt.axhline(0, c='k', ls='--', lw=0.5)
# plt.xlim(0, lk[-1])
# plt.ylim([-5, 5])
# plt.title(r'All bands (SZP)')
# plt.show()

E = np.linspace(-5, 5, 400)
bz = sisl.MonkhorstPack(H, [101, 1, 1])
plt.plot(E, bz.apply.average.DOS(E))
plt.xlim(E[0], E[-1]); plt.ylim(0, None)
plt.xlabel(r'$E - E_F$ [eV]')
plt.ylabel(r'DOS [1/eV]')
plt.legend()
plt.show()
