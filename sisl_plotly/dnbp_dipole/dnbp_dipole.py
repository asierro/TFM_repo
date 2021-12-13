import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize


def rodrigues(v, k, alpha):
    '''
    Return the rotated vector v about vector k by an angle alpha, using the Rodrigues formula
    '''
    k = k / np.linalg.norm(k)
    return v * np.cos(alpha) + np.cross(k, v) * np.sin(alpha) + k * k.dot(v) * (1 - np.cos(alpha))


def dipole_energy(dipole, Efield):
    """
    Calculates dipole energy
    :param dipole: array of vectors of dipoles
    :param Efield: electric field vector
    :return: array w/ -Efield*dipole for each dipole
    """
    return -np.apply_along_axis(np.dot, 1, dipole, Efield)


torsion = np.loadtxt('torsion_angles.txt')
fake_angle, energy = np.split(np.loadtxt('energies.dat'), 2, axis=1)
dipole = np.loadtxt('dipoles.dat')[:, 1:]
# megadipole = np.loadtxt('dipoles.dat')
symmetrize_angle = np.loadtxt('symmetrize_angle.txt')#[:, 2]

# Dipole vectors need to be rotated because the original molecule used was not aligned with the z axis (unequal z
# coordinate for outermost carbons), and later we are going to compare w/ molecules that ARE aligned
rot_vec = [-3.8046513925164684, 5.71929277405903, 0.0]
rot_ang = -0.32106480732727877

# # Save norm to check that it is conserved
# norm_dipole = np.apply_along_axis(np.linalg.norm, 1, dipole)

dipole = np.apply_along_axis(rodrigues, 1, dipole, rot_vec, rot_ang)
# new_norm_dipole = np.apply_along_axis(np.linalg.norm, 1, dipole)

# Rotate dipoles to be symmetrized and pointing up
sym_rot_vec = [6.0272876094560, 4.0095391340037, -0.0000000000000]
dipole_new = np.empty_like(dipole)
for i, dip in enumerate(dipole):
    for threesome in symmetrize_angle:
        if abs(threesome[0] - fake_angle[i]) < 0.1:
            dipole_new[i] = rodrigues(dip, sym_rot_vec, threesome[2]*np.pi/180.)


# Dipole atomic units to e*Ang: 1 a.u. = conversion_factor e*Ang
conversion_factor = 8.4783536255E-30 / 1.602E-29
dipole_new *= conversion_factor

# convert fake angles into real torsion angles using the torsion "lookup table"
real_angle = np.empty_like(fake_angle)
for i, fake in enumerate(fake_angle):
    for pair in torsion:
        if abs(pair[0]-fake) < 0.1:
            real_angle[i] = pair[1]

plus_angle = abs(real_angle)

# plt.scatter(plus_angle, energy); plt.show()

dip_en = dipole_energy(dipole_new, [0., 0., 0.25])
dip_en_minus = dipole_energy(dipole_new, [0., 0., -0.25])

# plt.scatter(plus_angle, energy.T+dip_en); plt.show()
# plt.scatter(plus_angle, energy.T+dip_en_minus); plt.show()

a=plt.scatter(list(plus_angle) + list(360 - plus_angle), list(energy.T + dip_en) + list(energy.T + dip_en_minus),
            marker='+')
# plt.show()

# Just the original 0 field energies:
b=plt.scatter(list(plus_angle) + list(360 - plus_angle), list(energy.T) + list(energy.T), marker='.')

# ----------------------------------------------------------------------------------------------------------------------
# Now let's plot the original scf calculation with 0.25eV/Ang electric field , to compare with this last plot

torsion2 = np.loadtxt('torsion_ang_0.25.txt')
fake_angle2, energy2 = np.split(np.loadtxt('scf_en_0.25.dat'), 2, axis=1)

# convert fake angles into real torsion angles using the torsion2 "lookup table"
real_angle2 = np.empty_like(fake_angle2)
for i, fake in enumerate(fake_angle2):
    for pair in torsion2:
        if abs(pair[0]-fake) < 0.1:
            real_angle2[i] = pair[1]

c=plt.scatter(real_angle2, energy2, marker='.')

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Palatino'], size=10)

ax = plt.gca()
ax.set_xlabel('Torsion angle (º)')
ax.set_ylabel(r'$E_{\mathrm{tot}}$ (eV)')
ax.legend([a,b,c],[r'Dipole approx. ($|\mathbf{E}|$=0.25V/$\textup{\AA}$)', r'SCF ($|\mathbf{E}|$=0)',
                   r'SCF ($|\mathbf{E}|$=0.25V/$\textup{\AA}$)'])

# plt.show()
plt.savefig('dnbp_dipole.pdf', format='pdf', dpi=300, bbox_inches='tight')


# ----------------------------------------------------------------------------------------------------------------------

# Find minimum energy angle (m.e.a.). Then find m.e.a vs E field

# For that we will use a quartic function to fit the data and find the minimum of that.

'''
def quartic(x, center_4, amp_4, center_2, amp_2, h):
    return amp_4 * (x - center_4)**4 + amp_2 * (x - center_2)**2 + h


# 0 field (nitros up)
# xx = plus_angle.flatten()
# yy = energy.flatten()
# # 0 field (nitros down) (its the same as nitros up)
# xx = list(360 - plus_angle)
# yy = list(energy.T)


def ang_min(Efield, up):
    dip_en = dipole_energy(dipole_new, [0., 0., Efield])
    dip_en_minus = dipole_energy(dipole_new, [0., 0., -Efield])

    if up:
        print('yawho')
        xx = plus_angle.flatten()
        yy = energy.flatten() + dip_en
    else:
        xx = 360. - plus_angle.flatten()
        yy = energy.flatten() + dip_en_minus

    y_sort = yy[xx.argsort()]
    x_sort = np.sort(xx)

    # Inital guess for params of the quartic f'n:
    c = np.min(y_sort)

    # A(x-a)^4 part:
    x1, y1 = x_sort[0], y_sort[0]
    x2, y2 = x_sort[-1], y_sort[-1]
    P4 = ((y2-c) / (y1-c))**(1/4)
    a = (P4 * x1 - x2)/(P4 - 1)
    A = (y1 - c) / (x1 - a)**4

    # B(x-b)^2 part:
    x3, y3 = x_sort[2], y_sort[2]
    x4, y4 = x_sort[-2], y_sort[-2]
    P2 = ((y4-c)/(y3-c))**(1/2)
    b = (P2 * x3 - x4) / (P2 - 1)
    B = (y3 - c) / (x3 - b)**2

    center_4_0 = a
    amp_4_0 = A / 2
    center_2_0 = b
    amp_2_0 = B / 2
    h_0 = c

    popt, pcov = curve_fit(quartic, x_sort, y_sort, p0=[center_4_0, amp_4_0, center_2_0, amp_2_0, h_0], maxfev=500000)

    # Minimum angle:
    fit = minimize(quartic, x0=np.array([95]), args=tuple(popt))
    min_ang = fit.x[0]

    # plt.plot(x_sort, quartic(x_sort, *popt))
    # plt.scatter(x_sort, y_sort)
    return min_ang

up_angmins = []
down_angmins = []
for efi in np.linspace(0, 0.75, 101):
    up_angmins.append(-ang_min(efi, True) + 90.1780721)
    down_angmins.append(-ang_min(efi, False) + 269.3383187)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Palatino'], size=16)

a, = plt.plot(np.linspace(0, 0.75, 101), up_angmins, linestyle='dashed', color='k')
b, = plt.plot(np.linspace(0, 0.75, 101), down_angmins, linestyle='dotted', color='k')
ax = plt.gca()

ax.set_ylabel(r'$-\Delta \alpha$ (º)')
ax.set_xlabel(r'$|\mathbf{E}|$ (V/$\textup{\AA}$)')

ax.set_xticks(np.linspace(0, 0.75, 6))

ax.legend([a,b], [r'$\alpha<180$º',r'$\alpha>180$º'])
plt.savefig('delta_ang.pdf', format='pdf', dpi=300, bbox_inches='tight')

plt.show()
'''

