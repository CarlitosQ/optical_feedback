'''
Created on Jun 5, 2018

@author: tqian
'''

import scipy
import numpy as np
from numpy import pi
import cmath
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

c = 3e8
WL = np.linspace(1548e-9, 1552e-9, 10000)
WL_D = 1550e-9
freq = c / WL
freq_center = c / (WL_D)
L_gain = 300e-6
L_wg = (180 + 170 + 250) * 1e-6
L_grat = 999.54 * 1e-6
L_grat_eff = L_grat / 2
L_poly = (180 + 170 + 250) * 1e-6 + L_grat_eff
# L_poly = (180 + 170 + 250) * 1e-6 + 4.34e-4
L_ext = 3000e-6
L_ext = 1550e-9 * 2 * 1e3
L_ext = 0
L_ext = np. linspace(0e-6, 8000e-6, 10000)
L_ext = 3559.965288e-6
L_ext = 6159.965288e-6
L_ext = 77.02951466e-3
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

L_ext_over_cav = (n_poly * L_ext) / (n_gain * L_gain + n_poly * L_poly)
print ("L_ext / L_cav = %s" % L_ext_over_cav)