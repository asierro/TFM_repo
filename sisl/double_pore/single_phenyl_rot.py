import matplotlib.pyplot as plt

with open('energies_00.txt', 'r') as file:
    data = file.readlines()

ang2_list = list(range(-75, 105, 15))
delta_E = []
for i in range(len(data)):
    Ecb1 = float(data[i].split()[1])
    Ecb2 = float(data[i].split()[2])
    delta_E.append(Ecb2-Ecb1)

plt.scatter(ang2_list, delta_E); plt.show()
