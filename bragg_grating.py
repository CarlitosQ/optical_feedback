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

c = 3e8
n_gain = 3.2
n_air = 1.
n_poly = 1.46
WL = np.linspace(1548e-9, 1552e-9, 10000)
WL_D = 1550e-9
freq = c / WL
freq_center = c / (WL_D)
L_grat = 700.644*1e-6

kappa_L = 0.68641
kappa = kappa_L / L_grat

delta_beta = 2 * pi * n_poly * (1 / WL - 1 / WL_D)
gama = np.sqrt(kappa**2 + (1j * delta_beta)**2)
r = 1j * kappa * np.sinh(gama * L_grat) / ((1j * delta_beta) *
                                           np.sinh(gama * L_grat) + gama * np.cosh(gama * L_grat))
r_bragg = np.absolute(r)
r_phi = (np.angle(r)) / pi

fig, ax1 = plt.subplots()
color = 'blue'
ax1.set_xlabel('Frequency [GHz]', fontsize=20)
ax1.set_ylabel('Reflectivity $R_{grating}$', fontsize=20)
x_axis = ((freq - freq_center) / 1e9)
y_axis = r_bragg ** 2
# ax1.plot(((freq - freq_center) / 1e9), r_bragg ** 2, label="$R_{bragg}$")
ax1.plot(x_axis, y_axis, label="$R_{grating}$", color='black')
ax1.tick_params(axis='both', which='both', top=False, right=False)
# ax1.grid()
ax1.plot(x_axis[6000], y_axis[6000], '#F14040', marker='o', markersize=10)
ax1.plot(x_axis[4000], y_axis[4000], '#515151', marker='o', markersize=10)
ax1.axvline(x=x_axis[6000], color='#F14040', linestyle='--')
ax1.axvline(x=x_axis[4000], color='#515151', linestyle='--')
fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.show()


# # kappa_dc=n*omega*ips*
# L = np.arange(0, 300, 0.05)
# kappa_dc = -31.91 - 75.33j
# kappa_ac = kappa_dc / 2
# delta = kappa_dc
# alpha = np.sqrt(np.square(abs(kappa_ac)) - np.square(delta))
# r = (-kappa_ac * np.sinh(alpha * L)) / \
#     (delta * np.sinh(alpha * L) - 1j * alpha * np.cosh(alpha * L))
# plt.figure(1)
# plt.plot(r)
# plt.show()
