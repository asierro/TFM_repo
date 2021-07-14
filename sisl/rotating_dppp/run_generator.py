with open('RUN.fdf', 'r') as file:
    data = file.readlines()

for ang_a in range(0, 95, 5):
    for ang_b in range(0, 95, 5):
        data[0] = 'SystemLabel rot_dppp_{}_{}\n'.format(ang_a, ang_b)
        data[8] = '%include STRUCT_{}_{}.fdf\n'.format(ang_a, ang_b)
        with open('runs/RUN_{}_{}.fdf'.format(ang_a, ang_b), 'w') as file:
            file.writelines(data)
