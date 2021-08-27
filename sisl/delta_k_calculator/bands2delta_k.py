import numpy as np
import matplotlib.pyplot as plt

band5_0_0 = np.loadtxt('band5_0_0.txt')[:, :2]
band6_0_0 = np.loadtxt('band6_0_0.txt')[:, :2]
band7_0_0 = np.loadtxt('band7_0_0.txt')[:, :2]
band8_0_0 = np.loadtxt('band8_0_0.txt')[:, :2]

band4_0_45 = np.loadtxt('band4_0_45.txt')[:, :2]
band5_0_45 = np.loadtxt('band5_0_45.txt')[:, :2]
band6_0_45 = np.loadtxt('band6_0_45.txt')[:, :2]
band7_0_45 = np.loadtxt('band7_0_45.txt')[:, :2]

band3_0_75 = np.loadtxt('band3_0_75.txt')[:, :2]
band4_0_75 = np.loadtxt('band4_0_75.txt')[:, :2]
band5_0_75 = np.loadtxt('band5_0_75.txt')[:, :2]
band6_0_75 = np.loadtxt('band6_0_75.txt')[:, :2]

band3_0_90 = np.loadtxt('band3_0_90.txt')[:, :2]
band4_0_90 = np.loadtxt('band4_0_90.txt')[:, :2]
band5_0_90 = np.loadtxt('band5_0_90.txt')[:, :2]
band6_0_90 = np.loadtxt('band6_0_90.txt')[:, :2]


bands_0_0 = np.stack((band5_0_0, band6_0_0, band7_0_0, band8_0_0), axis=0)
bands_0_45 = np.stack((band4_0_45, band5_0_45, band6_0_45, band7_0_45), axis=0)
bands_0_75 = np.stack((band3_0_75, band4_0_75, band5_0_75, band6_0_75), axis=0)
bands_0_90 = np.stack((band3_0_90, band4_0_90, band5_0_90, band6_0_90), axis=0)


# def continuous_band(n, k):
#     m = n - 5
#     return np.interp(k, bands[m, :, 0], bands[m, :, 1])
#
#
# # Generate bands with even energy steps
# step = 0.
# even_band = bands[0, 1, :]
# for i, en_i in enumerate(bands[0, :, 1]):
#


def inverse_band(n, en, bandsfile):
    '''
    Returns k value(s) of the band for a given energy
    '''
    matches = []
    diffs = []
    for i, en_i in enumerate(bandsfile[n, :, 1]):
        diff = np.abs(en_i - en)
        if diff < 0.01:  # This is the MAXIMUM tolerance (usually it will be lower)
            matches.append(bandsfile[n, i, 0])
            diffs.append(diff)
    if matches:  # if matches is not empty
        return matches[np.argmin(diffs)]
    return np.nan


def prepare_for_plot(vb_min, vb_max, cb_min, cb_max, bandsfile):
    en_arr1 = np.linspace(vb_min, vb_max, 501)
    en_arr2 = np.linspace(cb_min, cb_max, 501)
    k5_arr = np.empty_like(en_arr1)
    k6_arr = np.empty_like(en_arr1)
    k7_arr = np.empty_like(en_arr2)
    k8_arr = np.empty_like(en_arr2)

    for i in range(501):
        k5_arr[i] = inverse_band(0, en_arr1[i], bandsfile)
        k6_arr[i] = inverse_band(1, en_arr1[i], bandsfile)
        k7_arr[i] = inverse_band(2, en_arr2[i], bandsfile)
        k8_arr[i] = inverse_band(3, en_arr2[i], bandsfile)
    return (np.abs(np.concatenate((k6_arr - k5_arr, k7_arr - k8_arr))),
            np.concatenate((en_arr1, en_arr2)))


k_0_0, en_0_0 = prepare_for_plot(-1.21, -0.37, 0.37, 1.07, bands_0_0)
k_0_45, en_0_45 = prepare_for_plot(-1.39, -0.37, 0.37, 1.26, bands_0_45)
k_0_75, en_0_75 = prepare_for_plot(-1.46, -0.37, 0.37, 1.34, bands_0_75)
k_0_90, en_0_90 = prepare_for_plot(-1.45, -0.37, 0.37, 1.34, bands_0_90)
ax = plt.gca()
# ax.set_xlim(0., 0.025)
ax.set_aspect(aspect=0.05)
plt.scatter(k_0_0, en_0_0, s=1, marker='o')
plt.scatter(k_0_45, en_0_45, s=1, marker='o')
plt.scatter(k_0_75, en_0_75, s=1, marker='o')
#plt.scatter(k_0_90, en_0_90, s=1, marker='o')
plt.show()
