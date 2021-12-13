import sisl
import numpy as np
import matplotlib.pyplot as plt

'''
Creation of dp_pp structure with different twist angles of the phenyl groups (based off of the relaxed structure)
'''

XV_sile = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/dp_pp/relax/NPG_PP-nonplanar.XV')
dppp = sisl.Geometry.read(XV_sile)

# Move H atoms belonging to phenyls that were on the other side of the supercell
sc_y = dppp.cell[1]
dppp = dppp.move(-sc_y, atoms=19).move(-sc_y, atoms=24)
coords = dppp.axyz()

# Remove all H atoms
H_remove = []
for ia, a, idx_specie in dppp.iter_species():
    if idx_specie == 0:
        H_remove.append(ia)

dppp = dppp.remove(H_remove)
coords = dppp.axyz()

phenyl1L = [67, 68, 69, 70]
phenyl2L = [73, 74, 75, 76]
phenyl1R = [83, 82, 86, 80]
phenyl2R = [78, 85, 55, 89]

# rotated1L = dp_pp_ortho.rotate(61.33, coords[102] - coords[107], atoms=phenyl1L, origo=coords[107], only='xyz')
# rotated2L = dp_pp_ortho.rotate(118.23, coords[108] - coords[113], atoms=phenyl2L, origo=coords[113], only='xyz')
# rotated1R = dp_pp_ortho.rotate(61.33, coords[115] - coords[117], atoms=phenyl1R, origo=coords[117], only='xyz')
# rotated2R = dp_pp_ortho.rotate(118.23, coords[90] - coords[120], atoms=phenyl2R, origo=coords[120], only='xyz')

# sisl.plot(dppp, s=10, atom_indices=True); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# # ax.set_xlim(30., 40.)#; ax.set_ylim(0., 10.)
# # sisl.plot(rotated1R, s=10); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# plt.tight_layout(); plt.show()


def rotate_phenyl(geom, which, angle):
    '''
    Rotate a phenyl from the planar position, the rotation vectors go from right to left.
    geom must be based on dp_pp_ortho (same indices)
    '''
    assert which in ['1L', '2L', '1R', '2R']
    if which == '1L':
        rotated = geom.rotate(angle - 28.67, coords[66] - coords[71], atoms=phenyl1L, origo=coords[71], only='xyz')
    elif which == '2L':
        rotated = geom.rotate(angle + 28.23, coords[72] - coords[77], atoms=phenyl2L, origo=coords[77], only='xyz')
    elif which == '1R':
        rotated = geom.rotate(angle - 28.28, coords[79] - coords[81], atoms=phenyl1R, origo=coords[81], only='xyz')
    else:
        rotated = geom.rotate(angle + 28.85, coords[54] - coords[84], atoms=phenyl2R, origo=coords[84], only='xyz')
    return rotated


def rotate_dppp(geom, ang1L, ang2L, ang1R, ang2R):
    rp = rotate_phenyl(geom, '1L', ang1L)
    rp = rotate_phenyl(rp, '2L', ang2L)
    rp = rotate_phenyl(rp, '1R', ang1R)
    return rotate_phenyl(rp, '2R', ang2R)


# Create 1NN TB Hamiltonian (rotate beforehand, this is just for testing)
# MaxR is negative for some reason so replace the C Atoms in this geometry for C Atoms with the correct orbital radius
C = sisl.Atom(Z=6, R=1.6*1.01)
atoms = sisl.Atoms(atoms=[C], na=104)
dppp = sisl.Geometry(coords, atoms=atoms, sc=dppp.cell)
dppp.set_nsc([3, 3, 0])
ang_a = 0.; ang_b = 75.
rotating_dppp = rotate_dppp(dppp, ang_a, ang_a + ang_b, ang_a + 180, ang_a + ang_b + 180)

# We have the geometry. Now let's construct the TB model.

H_1nn = sisl.Hamiltonian(rotating_dppp)
H_1nn.construct([[0.1, 1.6], [0., -2.7]])  # Using construct doesn't take into account the orientation of pz orbitals

band = sisl.BandStructure(H_1nn, [[0.5, 0., 0.], [0., 0., 0.],[0., 0.5, 0.]], 301, ['X', r'$\Gamma$', 'Y'])
eigs = band.apply.array.eigh(gauge='r')

# Retrieve the tick-marks and the linear k points
xtick, xtick_label = band.lineartick()
lk = band.lineark()
plt.plot(lk, eigs)

plt.ylabel('Eigenspectrum [eV]')
plt.gca().set_ylim(-1., 3.)
plt.gca().xaxis.set_ticks(xtick)
plt.gca().set_xticklabels(xtick_label)

# Also plot x-major lines at the ticks
ymin, ymax = plt.gca().get_ylim()
for tick in xtick:
    plt.plot([tick, tick], [ymin, ymax], 'k')

plt.show()

# I will rotate one phenyl aº from the planar position from 0º to 90º,
# and the other bº from aº to 90º+aº.
# The other bridge will be rotated "antisymmetrically" in order for both ribbons to be equivalent.
# Specifically one ribbon is rotated onto the other with a 180º x-axis rotation

# for ang_a in range(0, 95, 5):
#     for ang_b in range(0, 95, 5):
#         rotating_dppp = rotate_dppp(dp_pp_ortho, ang_a, ang_a + ang_b, ang_a + 180, ang_a + ang_b + 180)
#         rotating_dppp.write('structures/STRUCT_{}_{}.fdf'.format(ang_a, ang_b))
        # f = plt.figure()
        # f.clear()
        # plt.close(f)
        # plt.clf()
        # sisl.plot(rotating_dppp, s=10); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
        # plt.tight_layout()
        # plt.savefig('/home/asier/Documents/master/tfm/sisl/rotating_dppp/gif/rot_{}_{}.png'.format(ang_a, ang_b))