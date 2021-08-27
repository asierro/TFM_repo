import numpy as np

band5_0_0 = np.loadtxt('band5_0_0.txt')
band6_0_0 = np.loadtxt('band6_0_0.txt')
band7_0_0 = np.loadtxt('band7_0_0.txt')
band8_0_0 = np.loadtxt('band8_0_0.txt')


# def f_band_0_0(n, k):
#     if n == 5:
#         return np.interp(k, band5_0_0[:, 0], band5_0_0[:, 1])
#     if n == 6:
#         return np.interp(k, band6_0_0[:, 0], band6_0_0[:, 1])
#     if n == 7:
#         return np.interp(k, band7_0_0[:, 0], band7_0_0[:, 1])
#     if n == 8:
#         return np.interp(k, band8_0_0[:, 0], band8_0_0[:, 1])
#     return np.empty(0)
#
#
# # delta_k = np.zeros_like(band5_0_0[:, :1])
# delta_k[:, 0] =
# max_delta_e = np.max(np.diff(band6_0_0[:, 1]))
# for i, k in enumerate(band5_0_0):
#     en5_i = band5_0_0[i, 1]
#     idx_tuple = np.where(np.abs(band6_0_0[:, 1]-en5_i) <= max_delta_e + 0.000001)
#     if idx_tuple[0].size != 0:

def make_invertible(band):
    # for decreasing (increasing: greater to lesser and minus to plus)
    copy = np.copy(band)
    counter = 0
    for i in range(500):
        if band[i + 1, 1] > band[i, 1]:
            copy[i, 1] = np.nan
            counter = 0
        elif np.abs(band[i + 1, 1] - band[i, 1]) < 0.000001:
            copy[i, 1] = band[i, 1] + 0.00000000001 * (50 - i)
            counter += 1
        else:
            copy[i, 1] = band[i, 1]
            counter = 0
    if band[500, 1] > band[499, 1]:
        copy[500, 1] = np.nan
    elif np.abs(band[500, 1] - band[499, 1]) < 0.000001:
        copy[500, 1] = copy[499, 1] - 0.00000000001 * 50
    indices = np.where(np.isnan(copy[:, 1]))
    return np.delete(copy, indices[0], 0)


inv_band5_0_0 = make_invertible(band5_0_0)
inv_band6_0_0 = make_invertible(band6_0_0)
inv_band7_0_0 = make_invertible(band7_0_0)
inv_band8_0_0 = make_invertible(band8_0_0)


def f_band_0_0(n, en):
    if n == 5:
        return np.interp(en, inv_band5_0_0[:, 1], inv_band5_0_0[:, 0])
    if n == 6:
        return np.interp(en, inv_band6_0_0[:, 1], inv_band6_0_0[:, 0])
    if n == 7:
        return np.interp(en, inv_band7_0_0[:, 1], inv_band7_0_0[:, 0])
    if n == 8:
        return np.interp(en, inv_band8_0_0[:, 1], inv_band8_0_0[:, 0])
    return np.empty(0)


kak = np.diff(inv_band5_0_0, axis=0)
