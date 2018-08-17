'''
Created on Nov 6, 2017

@author: tqian
'''

import scipy
import numpy as np
import matplotlib.pyplot as plt

#=========================================================================
# #=========================================================================
# # # equation suitable for strong feedback
# #=========================================================================
# n_g = 3.2
# n_pol = 1.45
# n_ext = 1.45
# L_g = 400
# # L_f2p: 180, L_p: 170, L_p2g: 250, L_grat/2: 999.54/2
# L_pol = 180 + 170 + 250 + 999.54 / 2
# L_ext = np.arange(1, 100, 0.1)      # unit in [mm]
# F_ext_dB = np.arange(-30, -5, 5)
# # F_ext_array = np.array([])
# # for f in F_ext_dB:
# #     F_ext = float(np.power(10, 0.1 * f))
# #     F_ext_array = np.append(F_ext_array, F_ext)
# # print F_ext_array
# # F_ext = [0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9]
# # c = 3e8
# # mu = c / 1550e-9
# # w = 2 * 3.14 * mu
# # w_tau = np.arange(0, 2 * 3.14, 0.01)
# fig = plt.figure()
# ax = plt.subplot(111)
# for ff in F_ext_dB:
#     f_ext = np.power(10, 0.1 * ff)
#     print f_ext
#     numerator = n_g * L_g + n_pol * L_pol
#     denominator = numerator + f_ext * n_ext * L_ext * 1000
#     linewidth_reduction_ratio = numerator / denominator * numerator / denominator
#     ax.plot(L_ext, linewidth_reduction_ratio,
#             label="%s dB" % ff)
# plt.xlabel('external cavity length [mm]')
# plt.ylabel(r'linewidth reduction ratio $\Delta \nu / \Delta\nu_0$')
# # Shrink current axis by 15%
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
# # Put a legend to the right of the current axis
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.show()
#=========================================================================

#=========================================================================
# # Another equation from Agrawal
#=========================================================================
n_g = 3.2
n_pol = 1.45
r_0 = (n_g - 1) / (n_g + 1)     # amplitude reflectivity
L_pol = 180 + 170 + 250 + 999.54 / 2
L_g = 400
L_ext = np.linspace(0, 100, 1000)  # unit in mm
c = 3e14    # unit in um
alpha = -2
tau = 2 * n_pol * L_ext * 1000 / c
tau_in = 2 * (n_g * L_g + n_pol * L_pol) / c

# f_ext = np.arange(0, 0.8, 0.1)
# for f in f_ext:
#     kappa = (1 - r_0 * r_0) * f / r_0
#     X = kappa * tau * np.sqrt(1 + alpha * alpha)
# #     F = (1 + alpha * alpha) / ((1 + X) * (1 + X))
#     F = (1 + X) * (1 + X)
#     plt.figure(1)
#     plt.plot(L_ext, F)

f_ext = 0.01
f_ext = 5e-5

kappa = (1 - r_0 * r_0) * np.sqrt(f_ext) / (r_0 * tau_in)
X = kappa * tau * np.sqrt(1 + np.square(alpha))
F = np.square(1 + X)

f_r = 8e9   # relaxation frequency 8GHz
L_fr = c / (2 * f_r) / 1000     # seperation of short and long cavity

plt.figure(1)
plt.plot(L_ext, F, label="f_ext %s" % f_ext)
plt.axvline(x=L_fr, linestyle='dashed')
plt.text(L_fr, 2, r'$\L_{f_r}$=%s' % L_fr, fontsize=15)
plt.title('Equation from Agrawal')
plt.xlabel('external cavity length [mm]')
plt.ylabel(r'max. reduction factor $F$')
plt.legend(loc=0)
# np.savetxt('test.out', F, delimiter=',')
# plt.show()

#=========================================================================
# # Equation from Petermann
#=========================================================================
n_g = 3.2
n_pol = 1.45
L_pol = 180 + 170 + 250 + 999.54 / 2
L_g = 400
L_ext = np.arange(0, 100, 0.1)  # unit in mm
alpha = -2
c = 3e14    # unit in um
# R = 0.2
R = np.square((n_g - 1) / (n_g + 1))      # square of amplitude reflectivity
C_e = (1 - R) / (2 * np.sqrt(R))
tau_ext = 2 * n_pol * L_ext * 1000 / c
tau_L = 2 * (n_g * L_g + n_pol * L_pol) / c
print ("tau_L=%s" % tau_L)
C = np.sqrt(f_ext) * 2 * C_e * tau_ext / tau_L * np.sqrt(1 + np.square(alpha))
factor = np.square(1 + C)

plt.figure(2)
plt.plot(L_ext, factor, label="f_ext %s" % f_ext)
plt.axvline(x=L_fr, linestyle='dashed')
plt.text(L_fr, 2, r'$\L_{f_r}$=%s' % L_fr, fontsize=15)
plt.title('Equation from Petermann')
plt.xlabel('external cavity length [mm]')
plt.ylabel('max. reduction factor')
plt.legend(loc=0)
plt.show()
