with open('RUN.fdf', 'r') as file:
    data = file.readlines()

for angle in [40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160]:
    data[0] = 'SystemLabel rot_dn_dppp_{0:03d}\n'.format(angle)
    data[6] = '%include STRUCT_{0:03d}.fdf\n'.format(angle)
    with open('runs/RUN_{0:03d}.fdf'.format(angle), 'w') as file:
        file.writelines(data)
