import numpy as np

#ar = np.random.randint(0,10,(6,6,3))
ar = np.empty((3,6,6), dtype=tuple)
for i in range(3):
    for j in range(6):
        for k in range(6):
            ar[i, j, k] = (i, j, k)

print(ar[2, 2, (0,1)])
