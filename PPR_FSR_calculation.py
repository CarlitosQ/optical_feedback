'''
Created on Mar 12, 2018

@author: charl
'''

import scipy
import numpy as np
from numpy import pi
import cmath
import math
import matplotlib.pyplot as plt
plt.style.use('classic')
plt.rcParams.update(
    {'font.family': 'Arial', 'font.size': 15, 'axes.linewidth': 2, 'xtick.major.width': 2, 'ytick.major.width': 2})
from matplotlib.ticker import FormatStrFormatter

c = 3e8
WL = np.linspace(1548e-9, 1552e-9, 10000)
WL_D = 1550e-9
freq = c / WL
freq_center = c / (WL_D)
L_gain = 300e-6
L_wg = (180 + 170 + 250) * 1e-6
# L_wg = (525) * 1e-6
# L_grat = 999.54 * 1e-6
L_grat = 700 * 1e-6
L_grat_eff = L_grat / 2
L_poly = L_wg + L_grat_eff
# L_poly = (180 + 170 + 250) * 1e-6 + 4.34e-4
# L_ext = 3000e-6
# L_ext = 1550e-9 * 2 * 1e3
# L_ext = 0
L_ext = np. linspace(0e-6, 8000e-6, 10000)
# L_ext = 3559.965288e-6
# L_ext = 6159.965288e-6
r1 = np.sqrt(0.9)
r_1 = 0.99
# r_ext = 0.183673469388
r_ext = 0.5
r_ext = 0
n_gain = 3.2
n_air = 1.
n_poly = 1.46
# n_clad = 1.45
# n_core = 1.48
n_eff = 1.46

kappa_L = 0.68641
kappa = kappa_L / L_grat

FSR_f = c / (2 * (n_gain * L_gain + n_poly * (L_poly + L_ext)))
print("FSR = %s GHz" % (FSR_f / 1e9))
# FSR_lambda = WL_D**2 / (2 * (n_gain * L_gain + n_poly * (L_poly + L_ext)))
# print("FSR = %s nm" % (FSR_lambda * 1e9))

plt.figure(1)
plt.xlabel('External cavity length [mm]', fontsize=20)
plt.ylabel('FSR [GHz]', fontsize=20)
x_axis = (L_ext * 1e3)
y_axis = (FSR_f / 1e9)
plt.plot(x_axis, y_axis)
marker_position = [3912, 4487, 4824, 5199, 6087, 6637, 7949]
for i in marker_position:
    plt.plot(x_axis[i], y_axis[i], '#F14040', marker='o', markersize=10)
    print(y_axis[i])
plt.show()
