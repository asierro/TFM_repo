import sisl
import sisl.viz
import numpy as np
import matplotlib.pyplot as plt

'''
Creation of dp_pp structure functionalized with two nitro groups, replacing the biphenyl bridges with
2,2'-dinitrobiphenyl, with different twist angles of the phenyl groups (while keeping nitro groups planar wrt the 
corresponding ring.
Based off of the relaxed planar structure.
'''

XV_sile = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/dn_dppp/DN_DPPP.XV')
dn_dppp = sisl.Geometry.read(XV_sile)


def plotter():
    sisl.plot(dn_dppp, s=10, atom_indices=True)
    ax = plt.gca(); ax.set_aspect('equal', adjustable='box')
    # ax.set_xlim(right=12); ax.set_ylim(top=10)
    plt.show()


plotter()
# Will continue in another environment where the lat version of sisl (with plotly, etc.) is installed