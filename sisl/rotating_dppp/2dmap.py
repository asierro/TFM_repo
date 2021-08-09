import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import scipy.stats

'''
2D map of Delta E at Gamma for both phenyl angles
'''

data = []
E_arr = np.zeros((19, 19))
for i, ang1 in enumerate(range(0, 95, 5)):
    with open('energies_{}.txt'.format(ang1), 'r') as file:
        data.append(file.readlines())
    for j in range(len(data[i])):
        Ecb1 = float(data[i][j].split()[1])
        Ecb2 = float(data[i][j].split()[2])
        E_arr[i, j] = Ecb2 - Ecb1

ang1_list = list(range(0, 95, 5))
ang2_list = list(range(0, 95, 5))

where_minx = np.argmin(E_arr, axis=1)
min_arrx = np.zeros_like(E_arr)
for i in range(19):
    min_arrx[i][where_minx[i]] = 1.

where_miny = np.argmin(E_arr, axis=0)
min_arry = np.zeros_like(E_arr)
for i in range(19):
    min_arry[where_miny[i]][i] = 1.

nzx = np.nonzero(min_arrx)
nzy = np.nonzero(min_arrx)
xx = list(nzx[0] * 5) + list(nzy[0] * 5)
yy = list(nzx[1] * 5) + list(nzy[1] * 5)

coeffs = np.polyfit(xx, yy, 1)
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xx, yy)


# plt.pcolormesh(ang2_list, ang1_list, E_arr, shading='gouraud', cmap='inferno')
# plt.pcolormesh(ang2_list, ang1_list, min_arrx + min_arry, shading='nearest', cmap='Greys')
# plt.scatter(yy, xx)
plt.plot([coeffs[1], coeffs[1] + 90 * coeffs[0]], [0, 90], color='k')
ax = plt.gca()
ax.set_xlim([0., 90.]); ax.set_ylim([0., 90.])
ax.set_aspect('equal', adjustable='box')

# x, y = np.meshgrid(ang2_list, ang1_list)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(x, y, E_arr, cmap=cm.coolwarm, antialiased=False)
plt.show()

# The way it is programmed rn, the x-axis of the plot is the angle of the 2nd phenyl (wrt 1st phenyl)
# (so y-axis is the angle of the 1st phenyl)
