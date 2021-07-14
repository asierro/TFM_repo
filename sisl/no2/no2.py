import sisl
import numpy as np
import matplotlib.pyplot as plt

'''
Creation of NO2 molecule for relaxation
'''

N = sisl.Atom(7)
O = sisl.Atom(8)

pos_list=[(25., 25., 25.), (25.9215, 24.6117, 25.), (24.0785, 24.6117, 25.)]

a_list = [N, O, O]
atoms = sisl.Atoms(atoms=a_list, na=3)
sc = sisl.SuperCell([50., 50., 50.,])
no2 = sisl.Geometry(pos_list, atoms=atoms, sc=sc)

# sisl.plot(no2); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# plt.tight_layout(); plt.show()

no2.write('STRUCT_NO2.fdf')