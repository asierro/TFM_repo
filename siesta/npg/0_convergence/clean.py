forcelist = []

with open('maxforce_iter2.txt', 'r') as forces:
    for line in forces:
        lsplit = line.split()
        if lsplit != []:
            forcelist.append(line.split()[1])

with open('cleanforces2.txt', 'w+') as clean:
    for force in forcelist:
        clean.write(str(force) + '\n')
        
