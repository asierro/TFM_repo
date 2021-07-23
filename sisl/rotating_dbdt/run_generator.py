with open('RUN.fdf', 'r') as file:
    data = file.readlines()

for angle in range(0, 180, 5):
    data[0] = 'SystemLabel rot_dbdt_{0:03d}\n'.format(angle)
    data[8] = '%include ZMATRIX_{0:03d}.fdf\n'.format(angle)
    with open('runs/RUN_{0:03d}.fdf'.format(angle), 'w') as file:
        file.writelines(data)
