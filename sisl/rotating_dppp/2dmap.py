import matplotlib.pyplot as plt
import numpy as np

'''
2D map of Delta E at Gamma for both phenyl angles
'''

data = []
E_arr = np.zeros((11, 19))
for i, ang1 in enumerate(range(0, 55, 5)):
    with open('energies_{}.txt'.format(ang1), 'r') as file:
        data.append(file.readlines())
    for j in range(len(data[i])):
        Ecb1 = float(data[i][j].split()[1])
        Ecb2 = float(data[i][j].split()[2])
        E_arr[i, j] = Ecb2 - Ecb1

ang1_list = list(range(0, 55, 5))
ang2_list = list(range(0, 95, 5))

plt.pcolormesh(ang2_list, ang1_list, E_arr, shading='nearest'); plt.show()

# The way it is programmed rn, the x-axis of the plot is the angle of the 2nd phenyl (wrt 1st phenyl)
# (so y-axis is the angle of the 1st phenyl)
