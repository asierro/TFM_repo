import numpy as np
import matplotlib.pyplot as plt

# THIS CALCULATES DELTA E AT X, NOT GAMMA
# torsion = np.loadtxt('torsion_dn_dppp.txt')
# vb_cb = np.loadtxt('vb_cb_en_AT_X.txt')
#
# a,=plt.plot(np.abs(torsion[:, 1] + torsion[:, 2]) / 2., vb_cb[:, 2] - vb_cb[:, 1])
# b,=plt.plot(np.abs(torsion[:, 1] + torsion[:, 2]) / 2., vb_cb[:, 4] - vb_cb[:, 3])
# ax = plt.gca()
# ax.legend([a,b], ['VB','CB'])
# plt.show()

torsion = np.loadtxt('torsion_dn_dppp.txt')
vb_en = np.loadtxt('vbs.txt')
cb_en = np.loadtxt('cbs.txt')

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Palatino'], size=16)

a, = plt.plot(np.abs(torsion[:, 1] + torsion[:, 2]) / 2., vb_en[:, 1] - vb_en[:, 0], linestyle='dashed', color='k',
              marker='o')
b, = plt.plot(np.abs(torsion[:, 1] + torsion[:, 2]) / 2., cb_en[:, 1] - cb_en[:, 0], linestyle='dotted', color='k',
              marker='o')
ax = plt.gca()

ax.set_ylabel(r'$\Delta E$ at $\Gamma$ (eV)')
ax.set_xlabel('Torsion angle (º)')

ax.legend([a,b], ['VB','CB'])
plt.savefig('delta_e.pdf', format='pdf', dpi=300, bbox_inches='tight')
# plt.show()
