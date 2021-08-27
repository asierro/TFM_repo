import sisl
import numpy as np
import matplotlib.pyplot as plt

coefs = np.polyfit([1.42, 2. * 1.42 * np.cos(30*np.pi/180), 2 * 1.42], [2.756, 0.071, 0.38], 3)
plt.plot(np.linspace(0, 10, 100), np.polyval(coefs, np.linspace(0, 10, 100)))
plt.show()