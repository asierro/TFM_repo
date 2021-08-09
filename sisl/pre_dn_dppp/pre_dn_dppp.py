import sisl
import numpy as np
import matplotlib.pyplot as plt

'''
Creation of DN-DPPP for later relaxation. Based off of the nonplanar DPPP structure, substituting the bridges with the
(previously relaxed) DNBP
Based off of the relaxed planar structure.
'''

XV_dppp = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/dp_pp/relax/NPG_PP-nonplanar.XV')
dp_pp = sisl.Geometry.read(XV_dppp)

coords = dp_pp.axyz()

bridge1 = [102, 103, 104, 107, 106, 105,
           16, 18, 19, 13,
           108, 109, 110, 113, 112, 111,
           0, 17, 14, 10]
bridge2 = [115, 119, 118, 117, 116, 122,
           22, 23, 26, 21,
           90, 114, 121, 120, 125, 91,
           25, 34, 27, 24]


gnr7_13 = dp_pp.remove(bridge1 + bridge2)
coords_7_13 = gnr7_13.axyz()

#------------------------------------------------------------------------------------------
# Onto the DNBP part

XV_dnbp = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/dbdt_65/dbdt.XV')
dnbp = sisl.Geometry.read(XV_dnbp)
dnbp = dnbp.remove([24, 25])  # Remove H atoms where the 7-13 will attach
coords_dnbp = dnbp.axyz()

vect = coords_dnbp[17] - coords_dnbp[16]  # "direction" of the DNBP

# First make the direction vector lie on the x-y plane
aux = np.copy(vect)
aux[2] = 0
ang = np.arccos(np.dot(vect, aux)/(np.linalg.norm(vect) * np.linalg.norm(aux))) * 180 / np.pi
dnbp = dnbp.rotate(ang, np.cross(vect, aux), origo=coords_dnbp[16], only='xyz')
coords_dnbp2 = dnbp.axyz()

# Rotate the whole thing 180º so that nitro groups are facing up. Not necessary but better for plots later
vect2 = coords_dnbp2[17] - coords_dnbp2[16]
dnbp = dnbp.rotate(180, vect2, origo=coords_dnbp2[16], only='xyz')
coords_dnbp3 = dnbp.axyz()

# Now, find out how much the first phenyl is rotated, and rotate the whole thing to set it up similarly
# to the original nonplanar relaxed DPPP.
c17x = coords_dnbp3[17][0]; c17y = coords_dnbp3[17][1]
c16x = coords_dnbp3[16][0]; c16y = coords_dnbp3[17][1]
c11x = coords_dnbp3[11][0]; c11y = coords_dnbp3[11][1]

dist = np.abs((c17x-c16x)*(c16y-c11y) - (c16x-c11x)*(c17y-c16y))/np.sqrt((c17x-c16x)**2 + (c17y-c16y)**2)

vect_perp_2d = np.array([0., 0.])
vect_perp_2d[0] = vect2[1]
vect_perp_2d[1] = -vect2[0]
vect_perp_2d = vect_perp_2d * dist / np.linalg.norm(vect_perp_2d)

point= np.array([0., 0., 0.])
point[0] = c11x - vect_perp_2d[0]
point[1] = c11y - vect_perp_2d[1]
point[2] = coords_dnbp3[16][2]

vect_ang = coords_dnbp3[11] - point
aux_ang = np.copy(vect_ang)
aux_ang[2] = 0
phenyl_ang = np.arccos(np.dot(vect_ang, aux_ang)/(np.linalg.norm(vect_ang) * np.linalg.norm(aux_ang))) * 180 / np.pi

# Rotate so that half a torsion angle is on one side and viceversa
rot_ang = 69.122/2 - phenyl_ang
dnbp = dnbp.rotate(-rot_ang, vect2, origo=coords_dnbp2[16], only='xyz')
dnbp = dnbp.move(-coords_dnbp3[16])
coords_dnbp4 = dnbp.axyz()

# Now append the DNBP in the bridge spots, in the middle of them with the correct bridge angle.
ang_dnbp = np.arccos(vect2[0]/np.linalg.norm(vect2))*180/np.pi
# Bridge 1:
bridge_vect1 = coords_7_13[77] - coords_7_13[53]
ang_bridge1 = np.arccos(bridge_vect1[0]/np.linalg.norm(bridge_vect1))*180/np.pi
dnbp1 = dnbp.rotate(ang_bridge1 - ang_dnbp, np.array([0., 0., 1.]), only='xyz')

len_bridge1 = np.linalg.norm(bridge_vect1)
len_dnbp = np.linalg.norm(vect2)
offset1 = coords_7_13[53] + bridge_vect1*(len_bridge1-len_dnbp)/2/len_bridge1
dn_dppp = gnr7_13.add(dnbp1, offset=offset1)

# Bridge 2:
# First we have to rotate the DNBP 180 degrees so that the ribbons are equivalent (see notebook)
dnbp2 = dnbp.rotate(180, vect2, only='xyz')

bridge_vect2 = coords_7_13[55] - coords_7_13[85]
ang_bridge2 = np.arccos(bridge_vect2[0]/np.linalg.norm(bridge_vect2))*180/np.pi
dnbp2 = dnbp2.rotate(-ang_bridge2 - ang_dnbp, np.array([0., 0., 1.]), only='xyz')

len_bridge2 = np.linalg.norm(bridge_vect2)
offset2 = coords_7_13[85] + bridge_vect2*(len_bridge2-len_dnbp)/2/len_bridge2
dn_dppp = dn_dppp.add(dnbp2, offset=offset2)

# dn_dppp = dn_dppp.tile(2, 0).tile(2, 1)
dn_dppp.write('STRUCT_DN_DPPP.fdf')

sisl.plot(dn_dppp, s=10); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# sisl.plot(gnr7_13, s=10, atom_indices=True); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# sisl.plot(dnbp2, atom_indices=True); ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
# ax.set_xlim(0., 10.); ax.set_ylim(-10., 10.)
plt.tight_layout()
plt.show()
