'''
Created on 13.11.2017

@author: tqian
'''

# import scipy
# import numpy as np
# import matplotlib.pyplot as plt

#=========================================================================
# #=========================================================================
# # small calculation
# #=========================================================================
# x = np.arange(0, 10, 0.1)
# alpha = -2
# F = (1 + alpha * alpha) / ((1 + x) * (1 + x))
# # F = ((1 + x) * (1 + x))
# plt.figure(1)
# plt.plot(x, F)
# plt.show()
#=========================================================================

#=========================================================================
# calculate A
#=========================================================================
# n_g = 3.2
# n_pol = 1.45
# c = 3e14      # unit in um
# v_g0 = c / n_g
# v_g1 = c / n_pol
# alpha = -2
# L_pol = 180 + 170 + 250 + 999.54 / 2
# L_g = 400
# L_0 = L_pol + L_g     # unit in um
# T_c = 0.5
# r_0 = 0.2
# A = v_g0 / (v_g1 * alpha * L_0) * (T_c / 2) / (r_0 + T_c / 2)
# print A

# def uniqueish_color():
#     """There're better ways to generate unique colors, but this isn't awful."""
#     return plt.cm.gist_ncar(np.random.random())
#
#
# xy = (np.random.random((10, 2)) - 0.5).cumsum(axis=0)
#
# print xy
# fig, ax = plt.subplots()
# for start, stop in zip(xy[:-1], xy[1:]):
#     x, y = zip(start, stop)
#     ax.plot(x, y, color=uniqueish_color())
# plt.show()


# for i in range(3, 0, -1):
#     print i

import sys

print (sys.version)