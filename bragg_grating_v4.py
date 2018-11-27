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
from matplotlib.ticker import FormatStrFormatter

c = 3e8
n_gain = 3.2
n_air = 1.
n_poly = 1.46
n_eff = 1.46
WL = np.linspace(1548e-9, 1552e-9, 10000)
# WL = np.linspace(1545e-9, 1555e-9, 1000000)
WL_D = 1550e-9
freq = c / WL
freq_center = c / (WL_D)
L_gain = 300e-6
# L_grat = 699.4075 * 1e-6
L_grat = 699.84 * 1e-6
L_grat = 700.644*1e-6
L_grat_eff = L_grat / 2
# 170 phase section L = 700e-6
L_poly = (180 + 190 + 155) * 1e-6 + L_grat_eff
L_poly = 509 * 1e-6 + L_grat_eff
# L_poly = (180 + 170 + 250) * 1e-6 + 4.34e-4
L_wg = 509 * 1e-6
# L_ext = 3209.75853784e-6
# L_ext = 6109.758538e-6
# L_ext = 1550e-9 * 2 * 1e3
# L_ext = 5060.19528784e-6
L_ext = (267 + 518.362787842 + 688) * 1e-6
L_ext = 6359.758538 * 1e-6
# 6359.758538
# 5309.758538
# 4869.758538

# print(L_ext + L_grat_eff)
# L_ext = 0

r1 = np.sqrt(0.9)
r_1 = 0.99
r_ext = (n_poly - n_air) / (n_poly + n_air)   # 0.18699186991869918
# r_ext = 0.5
# r_ext = 0

# n_clad = 1.45
# n_core = 1.48


kappa_L = 0.68641
kappa = kappa_L / L_grat

spacing = c / (2 * (n_gain * L_gain + n_poly * L_poly))
print("mode spacing for ADVA = %s GHz" % (spacing / 1e9))

external_cavity = c / (2 * n_poly * (L_ext + L_grat_eff))
print("external cavity spacing = %s GHz" % (external_cavity / 1e9))

FSR = c / (2 * (n_gain * L_gain + n_poly * (L_poly + (L_ext + L_grat_eff))))
print("FSR = %s GHz" % (FSR / 1e9))

gain = 15 * 100    # 27.368 cm^-1

# 1st way
# don't consider loss
delta_beta = 2 * pi * n_eff * (1 / WL - 1 / WL_D)
delta = delta_beta
gama = np.sqrt(kappa**2 + (1j * delta_beta)**2)

T11 = np.cosh(gama * L_grat) - 1j * delta / gama * np.sinh(gama * L_grat)
T12 = -1j * kappa / gama * np.sinh(gama * L_grat)
T21 = 1j * kappa / gama * np.sinh(gama * L_grat)
T22 = np.cosh(gama * L_grat) + 1j * delta / gama * np.sinh(gama * L_grat)
det_T = T11 * T22 - T21 * T12


r = 1j * kappa * np.sinh(gama * L_grat) / ((1j * delta_beta) *
                                           np.sinh(gama * L_grat) + gama * np.cosh(gama * L_grat))


r = kappa * np.sinh(gama * L_grat) / ((1j * delta_beta) *
                                      np.sinh(gama * L_grat) + gama * np.cosh(gama * L_grat))


# W1 = np.exp(-1j * (2 * pi * n_poly * L_wg / WL)) * np.exp(-8000 * L_wg)
W1 = np.exp(-2j * (2 * pi * n_poly * L_wg / WL))
# W1 = 1
W2 = np.exp(-2j * (2 * pi * n_poly * L_ext / WL))   # without propagation loss
W22 = np.exp(2j * (2 * pi * n_poly * L_ext / WL))
W3 = np.exp(-2j * (2 * pi * n_poly * L_wg / WL))
# gain = 27.368 * 100    # 27.368 cm^-1
beta_gain = 2 * pi * n_gain / WL + 1j * (gain) / 2
W_gain = np.exp(-2j * (beta_gain * L_gain))
W5 = W3 * W_gain
# r_eff = (r * W1 + r_ext * W2) / (1 + r * W1 * r_ext * W2)
r_eff = (r + r_ext * W2) / (1 + r * r_ext * W2) * W1
# r_eff = (r + r_ext * W2) / (1 + r * r_ext * W2)
# r_eff = r * W1
# r_eff = r + r_ext * (1 - r**2) * W2    # formula from Lang & Kobayashi

# =========================================================================
# r_eff = (np.cosh(gama * L_grat) - 1j * delta_beta / gama * np.sinh(gama *
#                                                                    L_grat) * r_ext * W1 * W2 + -1j * kappa / gama * np.sinh(gama * L_grat)) / (1j * kappa / gama * np.sinh(gama * L_grat) * r_ext * W1 * W2 + np.cosh(gama * L_grat) + 1j * delta_beta / gama * np.sinh(gama * L_grat))
# =========================================================================

T12 = -1j * kappa / gama * np.sinh(gama * L_grat)
T21 = 1j * kappa / gama * np.sinh(gama * L_grat)
r_eff = (T11 * r_ext * W1 * W2 + T12) / (T21 * r_ext * W1 * W2 + T22)

# r_eff = (T21 * W22 - r_ext * T22 * W2) / (T11 * W22 - r_ext * T12 * W2)
# # tqian

r_eff = T21 / T11 + (det_T / T11**2 * r_ext * W2) / \
    (1 + T12 / T11 * r_ext * W2)    # tqian derived

# =========================================================================
# r_eff = r + (det_T / T11**2 * r_ext * W2) / \
#     (1 + T12 / T11 * r_ext * W2)    # tqian derived
# =========================================================================


# gain = 27.368 * 100    # 27.368 cm^-1
# gain = 12 * 100    # 27.368 cm^-1
beta_3 = 2 * pi * n_gain / WL + 1j * (gain) / 2
G = r1 * np.exp(-2j * (beta_3 * L_gain)) * r_eff
G = r1 * np.exp(-2j * (beta_3 * L_gain)) * W1 * r_eff


# =========================================================================
# r_left = T21 / T11 + ((1 / T11 * det_T / T11 * r1 * W5) /
#                       (1 + T12 / T11 * r1 * W5))    # by tqian
#
# r_left = 1 / T11 * (T21 + (det_T * r1 * W5) /
#                     (T11 + T12 * r1 * W5))  # from literature
#
# r_right = r_ext * W2
# G = r_left * r_right
# =========================================================================

r_bragg = np.absolute(r)
r_phi = (np.angle(r)) / pi

G_abs = np.absolute(G)
# G_phi = G / G_abs
# G_phi = np.arccos(np.real(G) / G_abs) / pi
G_phi = np.angle(G) / (pi)

# L_eff = -(WL**2/(4*np.pi*n_poly))*(np.diff(r_phi)/np.diff(WL))
# print(L_eff)

# 1st plot
fig, ax1 = plt.subplots()
color = 'blue'
ax1.set_xlabel('Frequency (GHz)')
ax1.set_ylabel('Reflectivity $R_{bragg}$', color=color)
ax1.plot(((freq - freq_center) / 1e9), r_bragg **
         2, label="$R_{bragg}$", color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.grid()
# ax1.legend(loc=0)
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'green'
# we already handled the x-label with ax1
ax2.set_ylabel('Phase/$\pi$', color=color)
ax2.plot(((freq - freq_center) / 1e9), (r_phi),
         label="$\phi_{bragg}$", color=color)
ax2.tick_params(axis='y', labelcolor=color)
# ax2.legend(loc=0)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# plt.title("$r_{ext}$=%s, $L_{ext}$=%s mm" % (r_ext, L_ext * 1000))
fig.tight_layout()  # otherwise the right y-label is slightly clipped

# 2nd plot
freq_axis = ((freq - freq_center) / 1e9)
fig, ax1 = plt.subplots()
color = 'blue'
ax1.set_xlabel('Frequency (GHz)')
ax1.set_ylabel('Round Trip Gain $|G|$', color=color)
ax1.plot(((freq - freq_center) / 1e9), (G_abs**2), label="$|G|$", color=color)
ax1.set_xlim(-150, 150)
ax1.tick_params(axis='y', labelcolor=color)
ax1.grid()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'green'
# we already handled the x-label with ax1
ax2.set_ylabel('Phase/$\pi$', color=color)
ax2.plot(((freq - freq_center) / 1e9), (G_phi), label="$\phi$", color=color)
ax1.set_xlim(-150, 150)
ax2.tick_params(axis='y', labelcolor=color)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# ax2.grid()
plt.title("$r_{ext}$=%s, $L_{ext}$=%s mm" % (r_ext, L_ext * 1000))
fig.tight_layout()  # otherwise the right y-label is slightly clipped


# =========================================================================
# np.savetxt('RTG_r_ext_%s_%smm_L_ext.txt' % (r_ext, L_ext * 1e3), np.transpose(
#     [((freq - freq_center) / 1e9), G_abs, G_phi]), fmt='%1.8e', delimiter='\t')
# =========================================================================

plt.show()
