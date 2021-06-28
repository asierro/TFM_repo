import sisl
import numpy as np
import matplotlib.pyplot as plt

graphene = sisl.geom.graphene().tile(2, axis=0)
# sisl.plot(graphene); plt.show()
H = sisl.Hamiltonian(graphene)
H.construct([[0.1, 1.43], [0., -2.7]])
# H.write('ELEC.nc')

H_device = H.tile(3, axis=1)
# sisl.plot(H_device); plt.show()
# H_device.write('DEVICE.nc')
tbt = sisl.get_sile('siesta.TBT.nc')

# plt.plot(tbt.E, tbt.transmission(), label='k-averaged')
# plt.plot(tbt.E, tbt.transmission(kavg=tbt.kindex([0, 0, 0])), label=r'$\Gamma$')
# plt.xlabel('Energy [eV]'); plt.ylabel('Transmission'); plt.ylim([0, None]); plt.legend(); plt.show()

# plt.plot(tbt.E, tbt.DOS(), label='DOS')
# plt.plot(tbt.E, tbt.ADOS(), label='ADOS')
# plt.xlabel('Energy [eV]'); plt.ylabel('DOS [1/eV]'); plt.ylim([0, None]); plt.legend(); plt.show()

# for i, coord in enumerate(H_device.axyz()):
#     print(f'{i}, {coord}')

# Even atoms belong to a sublattice,
# odd atoms to the other

plt.plot(tbt.E, tbt.DOS(atoms=list(range(0, 12, 2))), label='DOS')
plt.plot(tbt.E, tbt.ADOS(atoms=list(range(0, 12, 2))), label='ADOS')
plt.xlabel('Energy [eV]'); plt.ylabel('DOS [1/eV]'); plt.ylim([0, None]); plt.legend(); plt.show()
# It's just half of the whole device DOS