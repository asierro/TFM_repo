import sisl
# import sisl.viz.plotly    # does not work
import numpy as np
import matplotlib.pyplot as plt

graphene = sisl.geom.graphene(orthogonal=True)
# sisl.plot(graphene); plt.show()

H = sisl.Hamiltonian(graphene)
H.construct([[0.1, 1.43], [0., -2.7]])

# # For that, use plotly (GeometryPlot) to see indices in animation
# H_device = H.tile(3, axis=1)
# # sisl.plot(H_device); plt.show()
#
# # H_device.plot()    # not working
# # for i, coord in enumerate(H_device.axyz()):
# #     print(f'{i}, {coord}')
#
# H_device[5, 5] = 0.2    # There's one orbital per atom
# H_device.write('DEVICE_c.nc')

tbt = sisl.get_sile('run_1000_c/siesta.TBT.nc')

plt.plot(tbt.E, tbt.transmission(), label='k-averaged')
plt.plot(tbt.E, tbt.transmission(kavg=tbt.kindex([0, 0, 0])), label=r'$\Gamma$')
plt.xlabel('Energy [eV]'); plt.ylabel('Transmission'); plt.ylim([0, None]); plt.legend(); plt.show()

# plt.plot(tbt.E, tbt.DOS(), label='DOS')
# plt.plot(tbt.E, tbt.ADOS(), label='ADOS')
# plt.xlabel('Energy [eV]'); plt.ylabel('DOS [1/eV]'); plt.ylim([0, None]); plt.legend(); plt.show()