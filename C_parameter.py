'''
Created on Nov 6, 2017

@author: tqian
'''

import scipy
import numpy as np
import matplotlib.pyplot as plt
from math import pi

n_g = 3.2
n_pol = 1.45
r_0 = (n_g - 1) / (n_g + 1)
c = 3e14    # unit in um
alpha = -2    # linewidth enhancement factor
L_ext = np.arange(1, 100, 0.1)  # unit in mm
L_pol = 180 + 170 + 250 + 999.54 / 2
L_g = 400
tau = 2 * L_ext * 1000 / c
tau_in = 2 * (n_g * L_g + n_pol * L_pol) / c
# print tau_in
r = (tau_in * r_0 * c) / ((1 - r_0 * r_0) *
                          2 * L_ext * 1000 * np.sqrt(1 + alpha * alpha))
r_dB = 10 * np.log10(r)
# C = 3 * pi / 2      # C parameter
# r1 = C * r
f_r = 8e9   # relaxation frequency 8GHz
L_fr = c / (2 * f_r) / 1000     # seperation of short and long cavity
plt.figure(1)
plt.plot(L_ext, r_dB)
# plt.plot(L_ext, r1)
plt.axvline(x=L_fr, linestyle='dashed')
plt.text(L_fr, -10, r'$\L_{f_r}$=%s' % L_fr, fontsize=15)
plt.xlabel('external cavity length [mm]')
plt.ylabel(r'feedback fraction $r$ [dB]')
plt.show()
