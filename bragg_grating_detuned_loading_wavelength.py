'''
Created on Mar 9, 2018

@author: tqian
'''

import scipy
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
plt.style.use('classic')
plt.rcParams.update(
    {'font.family': 'Arial', 'font.size': 15, 'axes.linewidth': 2, 'xtick.major.width': 2, 'ytick.major.width': 2})
from matplotlib.ticker import FormatStrFormatter

# ======================= Without feedback ===========================
c = 3e8
n_gain = 3.2
n_air = 1.
n_poly = 1.46
WL = np.linspace(1548e-9, 1552e-9, 10000)
WL_D = 1550e-9
freq = c / WL
freq_center = c / (WL_D)
L_gain = 300e-6
L_phase = 525 * 1e-6
L_grat = 700.644*1e-6
L_grat_eff = L_grat/2
L_poly = L_phase + L_grat_eff

# # L_ext = 0
# L_ext = (267 + 518.362787842 + 688) * 1e-6
# # L_ext = 6359.76*1e-6

# r_1 = 0.99

# r_ext = 0
# # r_ext = (n_poly - n_air) / (n_poly + n_air)   # 0.18699186991869918
# # r_ext = 0.6

# kappa_L = 0.68641
# kappa = kappa_L / L_grat
# print(kappa)


# delta_beta = 2 * pi * n_poly * (1 / WL - 1 / WL_D)
# gama = np.sqrt(kappa**2 + (1j * delta_beta)**2)
# r = 1j * kappa * np.sinh(gama * L_grat) / ((1j * delta_beta) *
#                                            np.sinh(gama * L_grat) + gama * np.cosh(gama * L_grat))
# r_bragg = np.absolute(r)
# r_phi = (np.angle(r)) / pi

# W1 = np.exp(-2j * (2 * pi * n_poly * L_phase / WL))
# W2 = np.exp(-2j * (2 * pi * n_poly * L_ext / WL))
# r_eff = (r + r_ext * W2) / (1 + r * r_ext * W2) * W1

# gain = 22 * 100    # 27.368 cm^-1
# # gain = 10 * 100    # 27.368 cm^-1

# beta_3 = 2 * pi * n_gain / WL + 1j * (gain) / 2
# G = r_1 * np.exp(-2j * (beta_3 * L_gain)) * r_eff
# r_bragg = np.absolute(r)
# r_phi = (np.angle(r)) / pi
# G_abs = np.absolute(G)
# G_phi = np.angle(G) / (pi)

# fig, ax1 = plt.subplots()
# color = 'blue'
# ax1.set_xlabel('Detuned Wavelength [nm]', fontsize=20)
# # ax1.set_ylabel('Reflectivity $R_{grating}$', fontsize=20)
# ax1.set_ylabel('Round Trip Gain $|G|$', fontsize=20)
# x_axis = ((WL - WL_D) * 1e9)
# # y_axis = r_bragg ** 2
# y_axis = G_abs ** 2
# # ax1.plot(((freq - freq_center) / 1e9), r_bragg ** 2, label="$R_{bragg}$")
# ax1.plot(x_axis, y_axis, label="$R_{grating}$", color='black')
# ax1.tick_params(axis='both', which='both', top=False, right=False)
# # ax1.grid()
# ax1.plot(x_axis[6000], y_axis[6000], '#F14040', marker='o', markersize=10)
# ax1.plot(x_axis[4000], y_axis[4000], '#515151', marker='o', markersize=10)
# ax1.axvline(x=x_axis[6000], color='#F14040', linestyle='--')
# ax1.axvline(x=x_axis[4000], color='#515151', linestyle='--')
# fig.tight_layout()  # otherwise the right y-label is slightly clipped

# plt.show()

# ======================= With feedback ===========================
c = 3e8
n_gain = 3.2
n_air = 1.
n_poly = 1.46
WL = np.linspace(1548e-9, 1552e-9, 10000)
WL_D = 1550e-9
freq = c / WL
freq_center = c / (WL_D)
L_gain = 300e-6
L_phase = 525 * 1e-6
L_grat = 700.644*1e-6
L_grat_eff = L_grat/2
L_poly = L_phase + L_grat_eff

L_ext = (267 + 518.362787842 + 688) * 1e-6

r_1 = 0.99

r_ext = 0
r_ext = (n_poly - n_air) / (n_poly + n_air)   # 0.18699186991869918

kappa_L = 0.68641
kappa = kappa_L / L_grat
# print(kappa)


delta_beta = 2 * pi * n_poly * (1 / WL - 1 / WL_D)
gama = np.sqrt(kappa**2 + (1j * delta_beta)**2)
r = 1j * kappa * np.sinh(gama * L_grat) / ((1j * delta_beta) *
                                           np.sinh(gama * L_grat) + gama * np.cosh(gama * L_grat))
r_bragg = np.absolute(r)
r_phi = (np.angle(r)) / pi

W1 = np.exp(-2j * (2 * pi * n_poly * L_phase / WL))
W2 = np.exp(-2j * (2 * pi * n_poly * L_ext / WL))
r_eff = (r + r_ext * W2) / (1 + r * r_ext * W2) * W1

gain = 10 * 100    # 27.368 cm^-1

beta_3 = 2 * pi * n_gain / WL + 1j * (gain) / 2
G = r_1 * np.exp(-2j * (beta_3 * L_gain)) * r_eff
r_bragg = np.absolute(r)
r_phi = (np.angle(r)) / pi
G_abs = np.absolute(G)
G_phi = np.angle(G) / (pi)

fig, ax1 = plt.subplots()
color = 'blue'
ax1.set_xlabel('Detuned Wavelength [nm]', fontsize=20)
ax1.set_ylabel('Round Trip Gain $|G|$', fontsize=20)
x_axis = ((WL - WL_D) * 1e9)
y_axis = G_abs**2
ax1.plot(x_axis, y_axis, label="$|G|$", color='black')
ax1.tick_params(axis='both', which='both', top=False, right=False)
ax1.plot(x_axis[5055], y_axis[5055], '#F14040', marker='o', markersize=10)
# ax1.plot(x_axis[4800], y_axis[4800], '#515151', marker='o', markersize=10)
ax1.axvline(x=x_axis[5055], color='#F14040', linestyle='--')
fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.show()
