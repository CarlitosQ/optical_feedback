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
# WL = np.linspace(1545e-9, 1555e-9, 1000000)
WL_D = 1550e-9
freq = c / WL
freq_center = c / (WL_D)
L_gain = 300e-6
L_grat = 999.54 * 1e-6
# L_grat = 6000 * 1e-6
L_grat_eff = L_grat / 2
L_poly = (180 + 170 + 250) * 1e-6 + L_grat_eff
L_wg = (180 + 170 + 250) * 1e-6
L_wg = 5000 * 1e-6
L_ext = 3000e-6
L_ext = 0
r1 = np.sqrt(0.9)
r_1 = 0.99
# r_ext = 0.183673469388
r_ext = 0
n_gain = 3.2
n_air = 1.
n_poly = 1.46
n_eff = 1.46

kappa_L = 0.68641
kappa = kappa_L / L_grat

kappa = 5 * 100
kappa_L = kappa * L_grat

spacing = c / (2 * (n_gain * L_gain + n_poly * L_poly))
print ("mode spacing for ADVA = %s GHz" % (spacing / 1e9))

passive_cavity = c / (2 * n_poly * (L_poly + L_ext))
print ("passive cavity spacing = %s GHz" % (passive_cavity / 1e9))

delta_f_z = c / (n_poly * L_grat * np.sqrt(1 + (kappa_L / pi)**2))
print ("delta_f_z = %s GHz" % (delta_f_z / 1e9))

delta_f_A = c / (2 * (n_gain * L_gain + n_poly * L_poly))
print ("delta_f_A = %s GHz" % (delta_f_A / 1e9))

delta_f_P = c / (2 * (n_poly * (L_poly + L_ext)))
print ("delta_f_P = %s GHz" % (delta_f_P / 1e9))

delta_f_T = c / (2 * (n_gain * L_gain + n_poly * (L_poly + L_ext)))
print ("delta_f_T = %s GHz" % (delta_f_T / 1e9))

R_BGA = delta_f_z / delta_f_A
print ("R_BGA = %s" % (R_BGA))

print ("delta_f_z / delta_f_P= %s" % (delta_f_z / delta_f_P))


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


W1 = np.exp(-2j * (2 * pi * n_poly * L_wg / WL))
W2 = np.exp(-2j * (2 * pi * n_poly * L_ext / WL))
r_eff = ((r + r_ext * W2) / (1 + r * r_ext * W2))

gain = 27.368 * 100    # 27.368 cm^-1
gain = 10 * 100
beta_3 = 2 * pi * n_gain / WL + 1j * (gain) / 2
G = r1 * np.exp(-2j * (beta_3 * L_gain)) * W1 * r
G = r1 * np.exp(-2j * (beta_3 * L_gain)) * W1 * r_eff


r_bragg = np.absolute(r)
r_phi = (np.angle(r)) / pi

G_abs = np.absolute(G)
G_phi = np.angle(G) / (pi)


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
# ax1.set_xlim(-50, 50)
ax1.tick_params(axis='y', labelcolor=color)
ax1.grid()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'green'
# we already handled the x-label with ax1
ax2.set_ylabel('Phase/$\pi$', color=color)
ax2.plot(((freq - freq_center) / 1e9), (G_phi), label="$\phi$", color=color)
# ax1.set_xlim(-50, 50)
ax2.tick_params(axis='y', labelcolor=color)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# ax2.grid()
plt.title("$r_{ext}$=%s, $L_{ext}$=%s mm" % (r_ext, L_ext * 1000))
fig.tight_layout()  # otherwise the right y-label is slightly clipped


#=========================================================================
# np.savetxt('RTG_r_ext_%s_%smm_L_ext.txt' % (r_ext, L_ext * 1e3), np.transpose(
#     [((freq - freq_center) / 1e9), G_abs, G_phi]), fmt='%1.8e', delimiter='\t')
#=========================================================================

plt.show()
