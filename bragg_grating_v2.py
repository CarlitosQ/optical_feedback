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

WL = np.linspace(1548e-9, 1552e-9, 10000)
# WL = 1555e-9
WL_D = 1550e-9
# print ("WL= %s" % WL)
L_gain = 300e-6
# L_gain = 120e-6
L_grat = 999.54 * 1e-6
L_poly = (180 + 170 + 250 + 999.54 / 2) * \
    1e-6    # 170 phase section L = 700e-6
# L_poly = 700e-6
L_ext = 400e-6
L_ext = 0
r1 = np.sqrt(0.99)
# r_ext = r * np.cos(4 * np.pi * L_ext / WL)
r_ext = np.sqrt(0.01)
r_ext = 0
n_gain = 3.2
n_air = 1.
n_poly = 1.46
# n_poly = 3.4
n_clad = 1.45
n_core = 1.48
n_eff = 1.46

#=========================================================================
# # n_0 = np.linspace(1.5, 1.6, 1000)
# # k = 2 * np.pi * n_0 / WL
# # # print ("k= %s" % k)
# # delta_k = k - np.pi / WL
# # delta_n = 0.1
# # a = 100e-6
# # V = 2 * np.pi / WL * a * np.sqrt(np.square(n_core) - np.square(n_clad))
# # yeta = 1 - 1 / np.square(V)
# # omega = np.pi * delta_n * yeta / WL
# # s = np.sqrt(np.square(omega) - np.square(delta_k))
# # # print ("s= %s" % s)
# # R = np.square(omega) * np.square(np.sinh(s * L)) / (np.square(delta_k)
# #                                                     * np.square(np.sinh(s * L)) + np.square(s) * np.square(np.cosh(s * L)))
#
# # delta = 2 * np.pi * n_eff * (1 / WL - 1 / WL_D)
# # delta_n_eff = 0.0004
# # theta = delta + 2 * np.pi / WL * delta_n_eff
#
# # alpha = np.sqrt(np.square(kappa) - np.square(theta) + 0j)
# # print ("alpha= %s" % alpha)
# # r = np.square(np.sinh(cons * L)) / \
# #     (np.square(np.cosh(cons * L) - np.square(theta) / np.square(kappa)))
#
# # r = -kappa * np.sinh(alpha * L) / (delta * np.sinh(alpha *
# # L) - 1j * alpha * np.cosh(alpha * L))
#
# # kappa = np.pi * delta_n_eff / WL
# # print ("kappa= %s" % kappa)
# # gama = np.sqrt(np.square(kappa) - np.square(delta) + 0j)
# # print ("gama= %s" % gama)
# # R = np.square(kappa * np.sinh(gama * L)) / (np.square(delta *
# # np.sinh(gama * L) + np.square(kappa * np.cosh(gama * L))))
#=========================================================================

c = 3e8
freq = c / WL
kappa = 1 / L_grat

# 1st way
# old
alpha = 0
delta_beta = 2 * pi * n_eff * (1 / WL - 1 / WL_D)
gama = np.sqrt(kappa**2 + (alpha + 1j * delta_beta)**2)
r = 1j * kappa * np.sinh(gama * L_grat) / ((alpha + 1j * delta_beta)
                                           * np.sinh(gama * L_grat) + gama * np.cosh(gama * L_grat))

# don't consider loss
delta_beta = 2 * pi * n_eff * (1 / WL - 1 / WL_D)
gama = np.sqrt(kappa**2 + (alpha + 1j * delta_beta)**2)
r = 1j * kappa * np.sinh(gama * L_grat) / ((-1j * delta_beta) *
                                           np.sinh(gama * L_grat) + gama * np.cosh(gama * L_grat))

#=========================================================================
# # 2nd way
# delta = 2 * pi * n_eff * (1 / WL - 1 / WL_D)
# delta_n_eff_avg = kappa * WL_D / pi
# ku = 2 * pi / WL_D * delta_n_eff_avg
# ku_bar = delta + ku
# gama_b = np.sqrt(kappa**2 - (1j * ku_bar)**2)
# r = 1j * (kappa / gama_b) * np.sinh(gama_b * L_grat) / (
#     np.cosh(gama_b * L_grat) - 1j * (ku_bar / gama_b) * np.sinh(gama_b * L_grat))
#=========================================================================

#=========================================================================
# # 3rd way (combine 1st and 2nd)
# delta = 2 * pi * n_eff * (1 / WL - 1 / WL_D)
# delta_n_eff_avg = kappa * WL_D / pi
# ku = 2 * pi / WL_D * delta_n_eff_avg
# ku_bar = delta + ku
# gama = np.sqrt(kappa**2 + (1j * ku_bar)**2)
# r = 1j * kappa * np.sinh(gama * L_grat) / (-(1j * ku_bar)
#                                            * np.sinh(gama * L_grat) + gama * np.cosh(gama * L_grat))
#=========================================================================

# plt.figure(1)
# plt.plot(WL, delta)
# plt.plot(WL, ku_bar)


#=========================================================================
# # 4th way
# delta = 2 * pi * n_eff * (1 / WL - 1 / WL_D)
# delta_n_eff_avg = kappa * WL / pi
# ku = 2 * pi / WL * delta_n_eff_avg
# ku_bar = delta + gama
# gama = np.sqrt(kappa**2 - ku_bar**2)
# r = np.square(np.sin(gama * L_grat)) / \
#     (np.square(np.cosh(gama * L_grat)) - np.square(ku_bar / kappa))
#=========================================================================

#=========================================================================
# # 5th way
# delta_beta = 2 * pi * n_eff * (1 / WL - 1 / WL_D)
# gama = np.sqrt(kappa**2 - delta_beta**2)
# R = kappa**2 * np.sinh(gama * L_grat) / (delta_beta**2 *
#                                          np.sinh(gama * L_grat) + kappa**2 * np.cosh(gama * L_grat))
#=========================================================================


W1 = np.exp(-1j * (2 * pi * n_poly * L_poly / WL))
W2 = np.exp(-1j * (2 * pi * n_poly * L_ext / WL))    # without propagation loss
r_eff = (r * W1 + r_ext * W2) / (1 + r * W1 * r_ext * W2)
# r_eff = (r + r_ext * W2) / (1 + r * r_ext * W2)
r_eff = r * W1
# r_eff = r + r_ext * (1 - r**2) * W    # formula from Lang & Kobayashi
G = r1 * np.exp(-1j * (2 * pi * n_gain * L_gain) / WL) * r_eff

# phase = 4 * np.pi * n_eff / (L * WL)

# r_phi = np.arctan(np.imag(r) / np.real(r))
# r_ref = np.real(r) / np.cos(r_phi)

r_bragg = np.absolute(r)
# r_phi = np.angle(r) * np.deg2rad
r_phi = np.arccos(np.real(r) / r_bragg) / pi
r_phi = np.arcsin(np.imag(r) / r_bragg) / pi
# r_phi = np.arctan2(np.imag(r), np.real(r)) / pi
r_phi = (np.angle(r)) / pi
# print (r_phi / r_phi2)

#=========================================================================
# # check phase calculation
# print (cmath.phase(complex(-1.0, -0.0)))
# print (np.angle(complex(-1.0, -0.0)))
# com = complex(-1.0, -0.0)
# print (math.atan2(np.imag(com), np.real(com)))
# print (np.arctan2(np.imag(com), np.real(com)))
#=========================================================================

r_eff_ref = np.absolute(r_eff)
r_eff_phi = r_eff / r_eff_ref
r_eff_phi = np.arccos(np.real(r_eff) / r_eff_ref) / pi
r_eff_phi = np.angle(r_eff) / pi

#=========================================================================
# value = ref_eff_phi / const_gc_phi
# print ("value=%s" % value)
#=========================================================================

G_abs = np.absolute(G)
G_phi = G / G_abs
G_phi = np.arccos(np.real(G) / G_abs) / pi
G_phi = np.angle(G) / pi

# spacing = c / (2 * (n_gain * L_gain + n_poly * L_poly))
# print ("mode spacing for ADVA = %s GHz" % (spacing / 1e9))

freq = c / (WL * 1e9)
freq_center = c / (WL_D * 1e9)

#=========================================================================
# plt.figure(1)
# # plt.plot((WL * 1e6), (np.real(r)), label="$\phi_{bragg}$")
# plt.plot((WL * 1e6), np.square(np.imag(r)), label="$r_{bragg}$")
# # plt.plot((WL * 1e6), np.square(np.imag(r_eff)), label="$r_{eff}$")
# plt.legend(loc=0)
# plt.xlabel('Wavelength [$\mu m$]')
# plt.ylabel('Reflectivity')
# # plt.plot(WL, np.square(np.imag(r)))
#
# # plt.figure(2)
# plt.plot((WL * 1e6), (np.real(r)), label="$\phi_{bragg}$")
# plt.legend(loc=0)
# # plt.plot(freq, np.square(np.imag(r)))
# # plt.plot(WL, np.imag(r))
# # plt.plot(WL, phase)
# plt.xlabel('Wavelength [$\mu m$]')
# plt.ylabel('Phase')
#=========================================================================


#=========================================================================
# plt.figure(1)
# plt.plot((WL * 1e9), r_bragg**2, label="$r_{bragg}$")
# # plt.plot((WL * 1e9), (r_phi), label="$\phi_{bragg}$")
# # plt.plot((WL * 1e9), (r_phi2), label="$\phi_{bragg}$")
# plt.plot((WL * 1e9), (r_phi4), label="$\phi_{bragg}$")
# plt.legend(loc=0)
#
#
# plt.figure(2)
# plt.plot((WL * 1e9), r_eff_ref**2, label="$r_{eff}$")
# plt.plot((WL * 1e9), r_eff_phi, label="$\phi_{eff}$")
# plt.legend(loc=0)
#
#
# plt.figure(3)
# plt.plot((WL * 1e9), G_abs, label="$|G|$")
# plt.plot((WL * 1e9), (G_phi), label="$\phi$")
# plt.legend(loc=0)
# plt.xlabel('Wavelength [$\mu m$]')
# plt.ylabel('Reflectivity')
#=========================================================================

# 1st plot
# (freq-freq_center)
fig, ax1 = plt.subplots()
color = 'blue'
ax1.set_xlabel('Frequency (GHz)')
ax1.set_ylabel('Reflectivity $R_{bragg}$', color=color)
ax1.plot((freq_center - freq), r_bragg**2, label="$R_{bragg}$", color=color)
ax1.tick_params(axis='y', labelcolor=color)
# ax1.legend(loc=0)
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'green'
# we already handled the x-label with ax1
ax2.set_ylabel('Phase/$\pi$', color=color)
ax2.plot((freq_center - freq), (r_phi), label="$\phi_{bragg}$", color=color)
ax2.tick_params(axis='y', labelcolor=color)
# ax2.legend(loc=0)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.title("$r_{ext}$=%s, $L_{ext}$=%s mm" % (r_ext, L_ext * 1000))
fig.tight_layout()  # otherwise the right y-label is slightly clipped


# 2nd plot
fig, ax1 = plt.subplots()
color = 'blue'
ax1.set_xlabel('Frequency (GHz)')
ax1.set_ylabel('Reflectivity $R_{eff}$', color=color)
ax1.plot((freq_center - freq), r_eff_ref**2, label="$R_{eff}$", color=color)
ax1.tick_params(axis='y', labelcolor=color)
# ax1.legend(loc=0)
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'green'
# we already handled the x-label with ax1
ax2.set_ylabel('Phase/$\pi$', color=color)
ax2.plot((freq_center - freq), r_eff_phi, label="$\phi_{eff}$", color=color)
ax2.tick_params(axis='y', labelcolor=color)
# ax2.legend(loc=0)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.title("$r_{ext}$=%s, $L_{ext}$=%s mm" % (r_ext, L_ext * 1000))
fig.tight_layout()  # otherwise the right y-label is slightly clipped


# 3rd plot
fig, ax1 = plt.subplots()
color = 'blue'
ax1.set_xlabel('Frequency (GHz)')
ax1.set_ylabel('Round Trip Gain $|G|$', color=color)
ax1.plot((freq_center - freq), G_abs, label="$|G|$", color=color)
ax1.tick_params(axis='y', labelcolor=color)
# ax1.legend(loc=0)
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'green'
# we already handled the x-label with ax1
ax2.set_ylabel('Phase/$\pi$', color=color)
ax2.plot((freq_center - freq), (G_phi), label="$\phi$", color=color)
ax2.tick_params(axis='y', labelcolor=color)
# ax2.legend(loc=0)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.title("$r_{ext}$=%s, $L_{ext}$=%s mm" % (r_ext, L_ext * 1000))
fig.tight_layout()  # otherwise the right y-label is slightly clipped

#=========================================================================
# # 4th plot
# fig, ax1 = plt.subplots()
# color = 'blue'
# ax1.set_xlabel('Frequency (GHz)')
# ax1.set_ylabel('Reflectivity $R_{bragg}$', color=color)
# ax1.plot((freq_center - freq), r_bragg**2, label="$R_{bragg}$", color=color)
# ax1.tick_params(axis='y', labelcolor=color)
# # ax1.legend(loc=0)
# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
# color = 'green'
# # we already handled the x-label with ax1
# ax2.set_ylabel('Reflectivity $R_{eff}$', color=color)
# ax2.plot((freq_center - freq), (r_eff_ref**2), label="$R_{eff}$", color=color)
# ax2.tick_params(axis='y', labelcolor=color)
# # ax2.legend(loc=0)
# ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#=========================================================================

#=========================================================================
# # 5th plot try
# fig, ax1 = plt.subplots()
# color = 'blue'
# ax1.set_xlabel('Frequency (GHz)')
# ax1.set_ylabel('Reflectivity $R_{bragg}$', color=color)
# ax1.plot((freq - freq_center), r_bragg**2, label="$R_{bragg}$", color=color)
# ax1.tick_params(axis='y', labelcolor=color)
# # ax1.legend(loc=0)
# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
# color = 'green'
# # we already handled the x-label with ax1
# ax2.set_ylabel('Reflectivity $R_{eff}$', color=color)
# ax2.plot((freq - freq_center), (r_eff_ref**2), label="$R_{eff}$", color=color)
# ax2.tick_params(axis='y', labelcolor=color)
# # ax2.legend(loc=0)
# ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#=========================================================================

# plt.title("$r_{ext}$=%s, $L_{ext}$=%s mm" % (r_ext, L_ext * 1000))
#=========================================================================
# plt.figure(3)
# plt.plot(L_ext, value)
#=========================================================================

#=========================================================================
# np.savetxt('RTG_%sum_gain_%s_%smm.txt' % (L_gain * 1e3, r_ext, L_ext * 1e3), np.transpose(
#     [(freq_center - freq), G_abs, G_phi]), fmt='%1.8e', delimiter='\t')
#=========================================================================

plt.show()
