with open('RUN.fdf', 'r') as file:
    data = file.readlines()

ang_a = 0
for j, ang_b in enumerate(range(-75, 105, 15)):
    data[0] = 'SystemLabel rot_dppp_{0:02d}_{1:02d}\n'.format(ang_a, j)
    data[8] = '%include STRUCT_{0:02d}_{1:02d}.fdf\n'.format(ang_a, j)
    with open('runs/RUN_{0:02d}_{1:02d}.fdf'.format(ang_a, j), 'w') as file:
        file.writelines(data)
