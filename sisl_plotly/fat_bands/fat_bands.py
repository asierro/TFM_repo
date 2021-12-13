import sisl
# import sisl.viz

XV_sile = sisl.get_sile('/home/asier/Documents/master/tfm/siesta/dndppp_xv/rot_dn_dppp_040.XV')
dn_dppp = sisl.Geometry.read(XV_sile)

hsx = sisl.io.siesta.hsxSileSiesta('/home/asier/Documents/master/tfm/siesta/dndppp_hsx/rot_dn_dppp_040.HSX')
HS = hsx.read_hamiltonian()
# H = sisl.physics.Hamiltonian.read(HSX_sile, geometry=dn_dppp)
