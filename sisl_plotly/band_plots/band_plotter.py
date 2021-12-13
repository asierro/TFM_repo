import sisl
import sisl.viz
import matplotlib.pyplot as plt

# bands = sisl.io.get_sile('/home/asier/Documents/master/tfm/siesta/dn_dppp_bands/DN_DPPP.bands')
for i in range(40, 170, 10):
    f = plt.figure()
    f.clear()
    plt.clf()
    plt.close(f)

    bands = sisl.io.get_sile(f'/home/asier/Documents/master/tfm/siesta/rot_dn_dppp/rot_dn_dppp_{i:03d}.bands')
    plot = bands.plot(backend='matplotlib')
    # print(plot)
    plot.update_settings(Erange=[-2., 3.])

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', serif=['Palatino'], size=12)

    locs, labels = plt.xticks()
    plt.xticks(locs, ['X', r'$\Gamma$', 'Y'])

    ax = plt.gca()
    ax.set_aspect(0.13)
    ax.xaxis.grid(True)
    ax.set_xlabel('')
    ax.set_ylabel('E (eV)')

    ax.legend([rf'$\sim${i}º'], loc='upper right')

    plt.savefig(f'dn_dppp_bands_{i:03d}.pdf', format='pdf', dpi=300)
# plot.show()
