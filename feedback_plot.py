'''
Created on 16.03.2018

@author: tqian
'''

import scipy
import numpy as np
import matplotlib.pyplot as plt

""" Calculation 1 """
WL = np.linspace(1540e-9, 1555e-9, 10000)
WL_D = 1550e-9
L = 1000e-6
n_eff = 1.45
cos = np.cos(4 * np.pi * L / WL)
theta = 2 * np.pi - 4 * np.pi * n_eff * (WL - WL_D) / np.square(WL_D)
theta = 4 * np.pi * n_eff * L / WL + 5 * np.sin(4 * np.pi * L / WL)
plt.figure(1)
plt.plot((WL * 1e6), theta)
# plt.plot(WL, cos)
plt.show()

#=========================================================================
# """ Calculation 2 """
# l = 0.0011
# R_1 = 0.3
# L = 100e-3
# n = 1.45
# r_a = np.square(n * l) * R_1 / np.square(L * (1 - R_1))
# r_b = np.square(3 * np.pi / 2) * r_a
# print ("r_a= %s" % (r_a))
# print ("r_b= %s" % (r_b))
#=========================================================================
