import matplotlib.pyplot as plt

with open('single_phenyl_rot.txt', 'r') as file:
    data = file.readlines()

angles = []
delta_E = []
for i in range(len(data)-1):
    angles.append(float(data[i].split()[0]))
    Ecb1 = float(data[i].split()[1])
    Ecb2 = float(data[i].split()[2])
    delta_E.append(Ecb2-Ecb1)

plt.scatter(angles, delta_E); plt.show()
