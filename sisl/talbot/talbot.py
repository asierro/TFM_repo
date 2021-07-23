import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv


def psi(n, y, k_c):
    return 1j**n * jv(n, 2 * k_c * y)


x_val = np.arange(-26, 26)
y_val = np.arange(0., 1200., 0.2)

x_arr, y_arr = np.meshgrid(x_val, y_val, indexing='ij')

width = 14 * 1.42 * np.cos(30 * np.pi / 180)

k = 0.008
prob_arr = np.abs(psi(x_arr, y_arr, k))**2

x_dist = x_arr / width  # x in angtrom instead of n units (1 n = width Ang)

plt.pcolormesh(x_dist, y_arr, prob_arr, shading='nearest'); plt.show()