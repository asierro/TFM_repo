import sisl
import numpy as np
import matplotlib.pyplot as plt

graphene = sisl.geom.graphene(orthogonal=True)

elec = graphene.tile(25, axis=0)
H = sisl.Hamiltonian(elec)
H.construct(([0.1, 1.43], [0., -2.7]))
H.write('ELEC.nc')

device = elec.tile(15, axis=1)
device = device.remove(
    device.close(
        device.center(what='cell'), R=10.)
)

# Each carbon atom in graphene has 3 nn, whichever atom has lt 3 is a dangling bond
dangling = [ia for ia in device.close(device.center(what='cell'), R=14.)
            if len(device.close(ia, R=1.43)) < 3]
# plt.gca().set_aspect('equal', adjustable='box')
# sisl.plot(device, s=1); plt.show()
device = device.remove(dangling)
# plt.gca().set_aspect('equal', adjustable='box')
# sisl.plot(device, s=1); plt.show()
edge = []
for ia in device.close(device.center(what='cell'), R=14.):
    if len(device.close(ia, R=1.43)) < 4:
        edge.append(ia)
edge = np.array(edge)

# Pretty-print the list of atoms (for use in sdata)
# Note we add 1 to get fortran indices
print(sisl.utils.list2str(edge + 1))

Hdev = sisl.Hamiltonian(device)
Hdev.construct(([0.1, 1.43], [0, -2.7]))
Hdev.geometry.write('device.xyz')
Hdev.write('DEVICE.nc')
