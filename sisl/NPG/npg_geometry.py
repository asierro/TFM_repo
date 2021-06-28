import sisl
import numpy as np
import matplotlib.pyplot as plt

graphene = sisl.geom.graphene(orthogonal=True)

# Create graphene sheet with the size of a NPG unit cell
sheet = graphene.tile(2, 0).tile(14, 1).rotate(-90, [0, 0, 1]).move([0, 7.1, 0])
# Align supercell correctly
sc_x = sheet.axyz()[111][0]
sc_y = -sheet.cell[0, 1]
sc_z = sheet.cell[2, 2]
sheet.set_supercell([sc_x, sc_y, sc_z])

# Remove extra carbon atoms
npg = sheet.remove([4, 7,    # left side
                    28, 29,  # left hole
                    77, 84,  # right hole
                    105, 106, 109, 110, 104, 107]) # right side

# Change dangling bonds to H atoms
coords = [[6.14878, 3.55, 0], [9.83805, 1.42, 0], [22.13561, 1.42, 0], [25.82488, 3.55, 0]]
H_list = []
for i in range(4):
    H_list += list(npg.close(coords[i], 1.43))
H = sisl.Atom(1)
npg.atoms.replace(H_list, H)

# npg.write('STRUCT.fdf')
#npg = npg.tile(2, 0).tile(2, 1)

plt.gca().set_aspect('equal', adjustable='box'); sisl.plot(npg, atom_indices=True); plt.show()
