import sisl
import numpy as np
import matplotlib.pyplot as plt


def torsion_angle(p0, p1, p2, p3):
    """
    Calculates torsion angle. Source:
    stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    """
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b1 /= np.linalg.norm(b1)

    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

# -> starting torsion angle is 91.16º

fdf = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/dbdt/RUN.fdf')
dbdt = sisl.Geometry.read(fdf, output=True)  # read from output files

rot_atoms = [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
all_atoms = list(range(28))

coords = dbdt.axyz()
# full_rot = dbdt.rotate(90, [1.0, 0.0, 0.0], atoms=all_atoms, origo=coords[0], only='xyz')

# # Rotate the dbdt from the relaxed position, the rotation vector goes from right to left.
# for angle in range(0, 185, 5):
#     rotated = dbdt.rotate(angle, coords[0] - coords[3], atoms=rot_atoms, origo=coords[3], only='xyz')
#     new_coords = rotated.axyz()
#     print(torsion_angle(new_coords[19], new_coords[14], new_coords[0], new_coords[5]))
#     rotated.write('structures/STRUCT_{}.fdf'.format(angle))

# sisl.plot(rotated, atom_indices=True); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# ax.set_xlim(right=20); ax.set_ylim(top=10)
# plt.show()

for k, angle in enumerate(range(0, 360, 5)):
    rotated = dbdt.rotate(angle, coords[0] - coords[3], atoms=rot_atoms, origo=coords[3], only='xyz')
    # f = plt.figure()
    # f.clear()
    # plt.close(f)
    # plt.clf()
    # sisl.plot(rotated); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
    # ax.set_xlim(right=20); ax.set_ylim(top=10)
    # plt.savefig('/home/asier/Documents/master/tfm/sisl/rotating_dbdt/gif/rot_{0:03d}.png'.format(angle))
    new_coords = rotated.axyz()
    print(torsion_angle(new_coords[19], new_coords[14], new_coords[0], new_coords[5]))
    rotated.write('structures/STRUCT_{0:03d}.fdf'.format(angle))
    with open('structures/STRUCT_{0:03d}.fdf'.format(angle), 'r') as struct_file, \
         open('zmatrices/ZMATRIX_{0:03d}.fdf'.format(angle), 'w') as zmat_file, \
         open('torsion_angles.txt') as torsion_file:

        struct_lines = struct_file.readlines()
        torsion_lines = torsion_file.readlines()

        zmat_lines = struct_lines[:8]
        zmat_lines.append('\n')
        zmat_lines.extend(struct_lines[40:])
        zmat_lines.append('\n')
        zmat_lines.extend(['ZM.UnitsLength Ang\n', 'ZM.UnitsAngle rad\n', '\n', '%block Zmatrix\n', 'molecule\n'])

        zmat_lines.append('1 0 0 0 ' + ' '.join(tuple(map(str, list(coords[19])))) + ' 1 1 1\n')  # 1st torsion atom

        v0 = coords[14] - coords[19]; norm_v0 = np.linalg.norm(v0)
        zmat_lines.append('1 1 0 0 ' + str(norm_v0) + ' ' + str(np.arctan2(np.linalg.norm(v0[:1]), v0[2])) + ' ' +
                          str(np.arctan2(v0[1], v0[0])) + ' 1 1 1\n')  # 2nd torsion atom

        v1 = -v0
        v2 = coords[0] - coords[14]; norm_v2 = np.linalg.norm(v2)
        zmat_lines.append('1 2 1 0 ' + str(norm_v2) + ' ' + str(np.arccos(np.dot(v1, v2)/(norm_v2 * norm_v0))) + ' 0.0 '
                          +
                          '1 1 1\n')  # 3rd torsion atom

        v3 = -v2
        v4 = coords[5] - coords[0]; norm_v4 = np.linalg.norm(v4)
        zmat_lines.append('1 3 2 1 ' + str(norm_v4) + ' ' + str(np.arccos(np.dot(v3, v4)/(norm_v4 * norm_v2))) + ' ANG ' +
                          '1 1 0\n')  # 4th atom. Constraint ('0') on the torsion angle 'ANG'

        zmat_lines.append('cartesian\n')
        for coord_line in range(10, 38):
            if coord_line not in [19 + 10, 14 + 10, 0 + 10, 5 + 10]:  # these atoms are in the molecule block
                Nspecies = struct_lines[coord_line][-10:-9]
                x_y_z = struct_lines[coord_line][1:-10]
                zmat_lines.append(Nspecies + ' ' + x_y_z + '1 1 1\n')

        zmat_lines.extend(['constants\n', 'ANG ' + str(float(torsion_lines[k]) * np.pi / 180) + '\n',
                           '%endblock Zmatrix\n'])

        zmat_file.writelines(zmat_lines)
