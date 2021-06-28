import sisl
import numpy as np
import matplotlib.pyplot as plt

# graphene = sisl.geom.graphene(orthogonal=True)
# # sisl.plot(graphene); plt.show()
#
# H = sisl.Hamiltonian(graphene)
# H.construct([[0.1, 1.43], [0., -2.7]])
# # print(H)
#
# H.write('ELEC.TSHS')

tbt = sisl.get_sile('run_500_TSHS/siesta.TBT.nc')

# plt.plot(tbt.E, tbt.transmission(), label='k-averaged')
# plt.plot(tbt.E, tbt.transmission(kavg=tbt.kindex([0, 0, 0])), label=r'$\Gamma$')
# plt.xlabel('Energy [eV]'); plt.ylabel('Transmission'); plt.ylim([0, None]); plt.legend(); plt.show()

plt.plot(tbt.E, tbt.DOS(), label='DOS')
plt.plot(tbt.E, tbt.ADOS(), label='ADOS')
plt.xlabel('Energy [eV]'); plt.ylabel('DOS [1/eV]'); plt.ylim([0, None]); plt.legend(); plt.show()
# Both methods yield the same result!