import numpy as np

with open('positions.txt', 'r') as mypresto_file:
    mypresto_lines = mypresto_file.readlines()
    for i in range(22):
        old_line = mypresto_lines[i].split()
        new_line = []
        vec = np.array([float(j) for j in old_line[2:5]])
        vec += 5 * np.ones_like(vec)  # move molecule up so it is w/in 1 sc
        vec = np.round(vec, 4)
        new_line.extend([str(k) for k in list(vec)])
        specie = old_line[1][0]
        if specie == 'C':
            new_line.append('1')
        elif specie == 'H':
            new_line.append('2')
        elif specie == 'N':
            new_line.append('3')
        elif specie == 'O':
            new_line.append('4')
        else:
            print('Error!!!')
        new_line.append('#')
        new_line.append(str(i+1) + ':')
        new_line.append(specie)
        print(' '.join(new_line))
