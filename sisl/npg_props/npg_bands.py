import sisl
import numpy as np
import matplotlib.pyplot as plt

# import os
# cwd = os.getcwd()

fdf = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/npg/4_bands/RUN.fdf')
H = fdf.read_hamiltonian()
print('Read!')

kpath = sisl.BandStructure(H, np.array([[0.5, 0., 0.], [0.] * 3, [0., 0.5, 0.]]), 100,
                           name=[r'$X$', r'$\Gamma$', r'$Y$'])
eig = kpath.apply.array.eigh().T  # calculate all eigenvalues

linear_k, k_tick, k_label = kpath.lineark(True)

Emin, Emax = -1, 4
plt.ylabel(r'$E-E_F$ [eV]')
plt.xlim(linear_k[0], linear_k[-1])
plt.xticks(k_tick, k_label)
plt.ylim(Emin, Emax)
for e in eig:
    plt.plot(linear_k, e, color='k')

plt.show()
