with open('RUN.fdf', 'r') as file:
    data = file.readlines()

data.extend(['\n', '%block ExternalElectricField\n', '0.000 0.000 0.250 V/Ang\n', '%endblock ExternalElectricField\n'])
# data.extend(['\n', 'MD.UseSaveZM T\n', 'DM.UseSaveDM T\n'])

# for angle in [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 135, 140, 145, 150, 155, 160, 165, 170,
#               175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 315, 320, 325, 330, 335,
#               340, 345, 350, 355]:
for angle in [30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]:
    data[0] = 'SystemLabel rot_dnbp_{0:03d}\n'.format(angle)
    data[8] = '%include ZMATRIX_{0:03d}.fdf\n'.format(angle)
    with open('runs_0.25/RUN_{0:03d}.fdf'.format(angle), 'w') as file:
    # with open('runs/RUN_{0:03d}.fdf'.format(angle), 'w') as file:
        file.writelines(data)
