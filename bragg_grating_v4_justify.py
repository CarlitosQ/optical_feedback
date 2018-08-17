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
WL = np.linspace(1548e-9, 1552e-9, 10000)
WL_D = 1550e-9
freq = c / WL
freq_center = c / (WL_D)
L_A = 140e-6
L_G = 780e-6
L_P = 0e-6
# L_ext = 0
r1 = np.sqrt(0.99)
r1 = 0.32
r_ext = np.sqrt(0.1)
# r_ext = 0.1
r_ext = 0
n_gain = 3.7
n_air = 1.
n_poly = 1.46
n_poly = 3.2
n_clad = 1.45
n_core = 1.48
n_eff = 3.244
n_eff = 3.7

# kappa_L = 1.8
# kappa = kappa_L / L_G

kappa = 20 * 100
kappa_L = kappa * L_G

delta_f_z = c / (n_gain * L_G * np.sqrt(1 + (kappa_L / pi)**2))
print ("delta_f_z = %s GHz" % (delta_f_z / 1e9))

delta_f_A = c / (2 * (n_gain * L_A + n_gain * L_G))
print ("delta_f_A = %s GHz" % (delta_f_A / 1e9))

#=========================================================================
# delta_f_P = c / (2 * (n_gain * (L_G + L_P)))
# print ("delta_f_P = %s GHz" % (delta_f_P / 1e9))
# print delta_f_z / delta_f_P
#=========================================================================

delta_f_T = c / (2 * (n_gain * L_A + n_eff * (L_G)))
print ("delta_f_T = %s GHz" % (delta_f_T / 1e9))

R_BGA = delta_f_z / delta_f_A
print ("R_BGA = %s" % (R_BGA))


# kappa = 30 * 100


# 1st way
delta_beta = 2 * pi * n_eff * (1 / WL - 1 / WL_D)
delta = delta_beta
gama = np.sqrt(kappa**2 + (1j * delta_beta)**2)

T11 = np.cosh(gama * L_G) - 1j * delta / gama * np.sinh(gama * L_G)
T12 = -1j * kappa / gama * np.sinh(gama * L_G)
T21 = 1j * kappa / gama * np.sinh(gama * L_G)
T22 = np.cosh(gama * L_G) + 1j * delta / gama * np.sinh(gama * L_G)
det_T = T11 * T22 - T21 * T12

r = 1j * kappa * np.sinh(gama * L_G) / ((1j * delta_beta) *
                                        np.sinh(gama * L_G) + gama * np.cosh(gama * L_G))

r = kappa * np.sinh(gama * L_G) / ((1j * delta_beta) *
                                   np.sinh(gama * L_G) + gama * np.cosh(gama * L_G))

beta_1 = 2 * pi * n_poly / WL
beta_2 = 2 * pi * n_poly / WL
W1 = np.exp(-2j * (beta_1 * L_G))
W2 = np.exp(-2j * (beta_2 * L_P))    # without propagation loss
# r_eff = (r * W1 + r_ext * W2) / (1 + r * W1 * r_ext * W2)

#=========================================================================
# r_eff = (r + r_ext * W2) / (1 + r * r_ext * W2) * W1
#
# r_eff = ((np.cosh(gama * L_G) - 1j * delta_beta / gama * np.sinh(gama * L_G)) * r_ext * W1 * W2 + (-1j * kappa / gama * np.sinh(gama * L_G))
#          ) / ((1j * kappa / gama * np.sinh(gama * L_G)) * r_ext * W1 * W2 + (np.cosh(gama * L_G) + 1j * delta_beta / gama * np.sinh(gama * L_G)))
#
# T12 = -1j * kappa / gama * np.sinh(gama * L_G)
# T21 = 1j * kappa / gama * np.sinh(gama * L_G)
# r_eff = (T11 * r_ext * W1 * W2 + T12) / (T21 * r_ext * W1 * W2 + T22)
#=========================================================================

r_eff = T21 / T11 + (det_T / T11**2 * r_ext * W2) / \
    (1 + T12 / T11 * r_ext * W2)    # tqian derived

#=========================================================================
# r_eff = r + (det_T / (T11**2) * r_ext * W2) / \
#     (1 + T12 / T11 * r_ext * W2)    # tqian derived
#=========================================================================

gain = 27.368 * 100     # 27.368 cm^-1
gain = 16 * 100
gain = 100 * 100
beta_3 = 2 * pi * n_gain / WL + 1j * (gain / 2)
G = r1 * np.exp(-2j * (beta_3 * L_A)) * r_eff
G = r1 * np.exp(-2j * (beta_3 * L_A)) * r


r_bragg = np.absolute(r)
r_phi = (np.angle(r)) / pi

r_eff_ref = np.absolute(r_eff)
r_eff_phi = np.angle(r_eff) / pi

G_abs = np.absolute(G)
G_phi = np.angle(G) / pi

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
plt.title("$r_{ext}$=%s, $L_{ext}$=%s mm" % (r_ext, L_P * 1000))
fig.tight_layout()  # otherwise the right y-label is slightly clipped

#=========================================================================
# # 2nd plot
# fig, ax1 = plt.subplots()
# color = 'blue'
# ax1.set_xlabel('Frequency (GHz)')
# ax1.set_ylabel('Reflectivity $R_{eff}$', color=color)
# ax1.plot(((freq - freq_center) / 1e9), r_eff_ref **
#          2, label="$R_{eff}$", color=color)
# # ax1.set_xlim(-100, 100)
# ax1.tick_params(axis='y', labelcolor=color)
# # ax1.legend(loc=0)
# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
# color = 'green'
# # we already handled the x-label with ax1
# ax2.set_ylabel('Phase/$\pi$', color=color)
# ax2.plot(((freq - freq_center) / 1e9), r_eff_phi,
#          label="$\phi_{eff}$", color=color)
# # ax2.set_xlim(-100, 100)
# ax2.tick_params(axis='y', labelcolor=color)
# # ax2.legend(loc=0)
# # ax2.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
# plt.title("$r_{ext}$=%s, $L_{ext}$=%s mm" % (r_ext, L_P * 1000))
# fig.tight_layout()  # otherwise the right y-label is slightly clipped
#=========================================================================


# 3rd plot
(freq - freq_center)
fig, ax1 = plt.subplots()
color = 'blue'
ax1.set_xlabel('Frequency (GHz)')
ax1.set_ylabel('Round Trip Gain $|G|$', color=color)
ax1.plot(((freq - freq_center) / 1e9), (G_abs**2), label="$|G|$", color=color)
ax1.set_xlim(-100, 100)
ax1.tick_params(axis='y', labelcolor=color)
ax1.grid()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'green'
# we already handled the x-label with ax1
ax2.set_ylabel('Phase/$\pi$', color=color)
ax2.plot(((freq - freq_center) / 1e9), (G_phi), label="$\phi$", color=color)
ax1.set_xlim(-100, 100)
ax2.tick_params(axis='y', labelcolor=color)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.title("$r_{ext}$=%s, $L_{ext}$=%s mm" % (r_ext, L_P * 1000))
fig.tight_layout()  # otherwise the right y-label is slightly clipped


path = "//hhi.de/benutzer/home/tqian/Master_thesis_tqian/calculations/r_bragg_RTG_verify/"
np.savetxt(path + 'RTG_squared_%sum_gain_%s_%smm_verify.txt' % (L_A * 1e3, r_ext, L_P * 1e3), np.transpose(
    [((freq - freq_center) / 1e9), G_abs**2, G_phi]), fmt='%1.8e', delimiter='\t')

plt.show()
