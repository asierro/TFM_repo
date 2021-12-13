import sisl
import matplotlib.pyplot as plt
from rotating_dnbp import torsion_root

def plotter():
    sisl.plot(dnbp, atom_indices=True); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(right=12); ax.set_ylim(top=10)
    plt.show()


fdf = sisl.get_sile('STRUCT_aligned.fdf')
dnbp = sisl.Geometry.read(fdf, output=True)  # read from output files
coords = dnbp.axyz()
# plotter()

with open('STRUCT_aligned.fdf', 'r') as struct_file, \
        open('ZMATRIX_aligned.fdf', 'w') as zmat_file:
    struct_lines = struct_file.readlines()

    zmat_lines = struct_lines[:8]
    zmat_lines.append('\n')
    zmat_lines.extend(struct_lines[38:])
    zmat_lines.append('\n')
    zmat_lines.extend(['ZM.UnitsLength Ang\n', 'ZM.UnitsAngle rad\n', '\n', '%block Zmatrix\n'])

    zmat_lines.append('cartesian\n')
    for coord_line in range(10, 36):
        Nspecies = struct_lines[coord_line][-10:-9]
        x_y_z = struct_lines[coord_line][1:-10]
        if coord_line in [26, 27]:
            zmat_lines.append(Nspecies + ' ' + x_y_z + '1 1 0\n')
        else:
            zmat_lines.append(Nspecies + ' ' + x_y_z + '1 1 1\n')

    zmat_lines.append('%endblock Zmatrix\n')

    zmat_file.writelines(zmat_lines)
