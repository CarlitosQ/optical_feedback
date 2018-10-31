'''
Created on Dec 15, 2017

@author: tqian
'''

import scipy
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import timeit

start = timeit.default_timer()

# =========================================================================
# # Define basic paramters
# =========================================================================
c = 3e8
n_gain = 3.2
n_poly = 1.46
n_air = 1.
L_poly = (180 + 170 + 250 + 999.54 / 2) * 1e-6    # 170 phase section
L_poly = (550 + 699.54 / 2) * 1e-6    # 170 phase section
L_gain = 500e-6
L_ext = np.linspace(0, 170e-3, 100000)
# L_ext = (267 + 518.362787842 + 688) * 1e-6
# L_ext = 3
alpha = -2
WL = np.linspace(1540e-9, 1570e-9, 10000)

# =========================================================================
# Define r_ext
# =========================================================================
# r_ext = 0.0078        # weak
# r_ext = 0.8        # very strong
r_ext = (n_poly - n_air) / (n_poly + n_air)   # 0.183673469388
# r_ext = 0.9
R_ext = np.square(r_ext)
print("R_ext = %s" % float('%.4g' % R_ext))

R_bragg = 0.3
r_bragg = np.sqrt(R_bragg)
tau_ext = 2 * n_poly * L_ext / c
tau_L = 2 * (n_gain * L_gain + n_poly * L_poly) / c
# print (tau_ext / tau_L)
# print ("tau_L= %s" % float('%.4g' % tau_c))

ratio = R_ext / R_bragg
print("R_ext/R_bragg= %s" % float('%.4g' % ratio))

kappa_ext = float(r_ext / r_bragg * (1 - np.square(r_bragg)))
print("kappa_ext= %s" % float('%.4g' % kappa_ext))

kappa = (1 - r_bragg) / tau_L * np.sqrt(r_ext / r_bragg)
print("kappa= %s" % ('%.6e' % Decimal(str(kappa))))

X = tau_ext / tau_L * kappa_ext
# =========================================================================
# plt.figure(2)
# plt.plot(L_ext * 1000, X, 'r', label="$\kappa_{ext}$ %s, $R_{ext}$ %s" %
#          (round(kappa_ext, 3), round(R_ext, 5)))
# =========================================================================
# C = X * np.sqrt(1 + np.square(alpha))
# print ("C= %s" % C)

f_r = 9e9   # relaxation frequency 9 GHz
# seperation of short and long cavity, in mm
L_fr = c / (2 * n_poly * f_r) * 1000
print("L_fr= %s mm" % L_fr)

""" 1. Weak feedback condition """
""" 1.1 r of pol and air interface """
# F_weak = 1 + C
# F_weak_sqr= np.square(F_weak)

# plt.figure(1)
# plt.plot(L_ext, F_weak_sqr, label="f_ext %s, weak condition" % r2_ext)
# plt.title('$F^2$ factor for weak feedback')
# plt.xlabel('external cavity length [mm]')
# plt.ylabel('max. reduction factor $F^2$')
# plt.legend(loc=0)

# plt.figure(2)
# plt.plot(L_ext, F_weak, label="f_ext %s" % r2_ext)
# plt.title('$F$ factor for weak feedback')
# plt.xlabel('external cavity length [mm]')
# plt.ylabel('max. reduction factor $F$')
# plt.legend(loc=0)
# plt.show()

""" 1.2 Different r2_ext """
# r2_ext = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]
#
# for r in r2_ext:
#     kappa_ext = r / r2s * (1 - np.square(np.abs(r2s)))
#     X = tau_ext / tau_L * kappa_ext
#     C = X * np.sqrt(1 + np.square(alpha))
#     F_weak_sqr = np.square(1 + C)
#     plt.figure(1)
#     plt.plot(L_ext, F_weak_sqr, label="f_ext %s" % r)
# plt.title('$F^2$ factor for weak feedback')
# plt.xlabel('external cavity length [mm]')
# plt.ylabel('reduction factor $F^2$')
# plt.legend(loc=0)
# plt.show()

""" 1.3 Different r2s """
# r2s = [0.1, 0.2, 0.3]
#
# for r in r2s:
#     kappa_ext = r2_ext / r * (1 - np.square(np.abs(r)))
#     X = tau_ext / tau_L * kappa_ext
#     C = X * np.sqrt(1 + np.square(alpha))
#     F_weak_sqr = np.square(1 + C)
#     plt.figure(1)
#     plt.plot(L_ext, F_weak_sqr, label="r2s= %s" % r)
# plt.title('$F^2$, weak feedback, different r2s')
# plt.xlabel('external cavity length [mm]')
# plt.ylabel('reduction factor $F^2$')
# plt.legend(loc=0)
# plt.show()

""" 2. Strong feedback condition """
""" 2.1 r of pol and air interface """
# F_ext = R2_ext
# print ("F_ext= %s" % F_ext)
# A = tau_ext / tau_L
# F_strong = 1 + F_ext * A
# # F_strong = 1 + X
# F_strong_sqr = np.square(F_strong)
#
# plt.figure(1)
# plt.plot(L_ext, F_strong_sqr, label="kappa= %s, strong condition" % kappa_ext)
# # plt.title('$F^2$, strong & weak formula')
# plt.title('$F^2$ strong feedback condition')
# plt.xlabel('external cavity length [mm]')
# plt.ylabel('max. reduction factor $F^2$')
# plt.legend(loc=0)
# plt.figure(2)
# plt.plot(L_ext, F_strong, label="f_ext %s" % f_ext)
# plt.title('$F$ factor for strong feedback')
# plt.xlabel('external cavity length [mm]')
# plt.ylabel('max. reduction factor $F$')
# plt.legend(loc=0)
# plt.show()

""" 2.2 Different f_ext (r) """
# f_ext = [0.1, 0.2, 0.25, 0.3, 0.4, 0.5]
#
# for f in f_ext:
#     A = tau_ext / tau_L
#     F_strong_sqr = np.square(1 + np.sqrt(f) * A)
#     plt.figure(2)
#     plt.plot(L_ext, F_strong_sqr, label="f_ext %s" % f)
# plt.title('$F^2$ factor for strong feedback')
# plt.xlabel('external cavity length [mm]')
# plt.ylabel('reduction factor $F^2$')
# plt.legend(loc=0)
# plt.show()

""" 3 Combined plot """
F_factor = []
# =========================================================================
# Define C_limit: for single mode, C is smaller than 1 or 3pi/2
# =========================================================================
C_limit = 1
# C_limit = 0.977 * 3 * np.pi / 2
# print("C_limit= %s" % C_limit)

for l_ext in L_ext:
    tau_EXT = 2 * n_poly * l_ext / c
    X = tau_EXT / tau_L * kappa_ext
    C = X * np.sqrt(1 + np.square(alpha))
#     print ("C= %s, l_ext= %s" % (C, l_ext))
    if C <= C_limit:
        L_c = l_ext * 1000
#         print ("L_c= %s mm" % L_c)
    if C > C_limit or kappa_ext > 0.1 or R_ext > 0.1:    # strong condition
        yeta_couple = 0.5
        F_ext = R_ext * yeta_couple
#         print ("strong condition F_ext=%s" % F_ext)
        A = tau_EXT / tau_L
        F = 1 + F_ext * A
        F_sqr = np.square(F)
        F_factor.append(F_sqr)
    elif C <= C_limit and kappa_ext < 0.1:       # weak condition
        #         F = 1 + C * np.cos(2 * np.pi * c / WL * tau_ext + np.arctan(alpha))
        F = 1 + C
        F_sqr = np.square(F)
        F_factor.append(F_sqr)
    else:
        F_sqr = 0
        F_factor.append(F_sqr)
# print (zip(L_ext, F_factor))
# print len(F_factor)
plt.figure(1)
if kappa_ext > 0.1 or R_ext > 0.1:
    plt.plot(L_ext * 1000, F_factor, 'r', label="$\kappa_{ext}$ %s, $R_{ext}$ %s" %
             (round(kappa_ext, 3), round(R_ext, 5)))
    plt.title('Strong condition, $F^2$ factor, mirror / polymer interface')
#     np.savetxt('F_factor_Strong_kappa%s.txt' % round(kappa_ext, 3), np.transpose(
#         [(L_ext * 1000), F_factor]), fmt='%1.8e', delimiter='\t')
elif kappa_ext < 0.01:
    plt.plot(L_ext * 1000, F_factor, 'b', label="$\kappa_{ext}$ %s, $R_{ext}$ %s" %
             (round(kappa_ext, 3), round(R_ext, 5)))
#     plt.title('Weak condition, $F^2$ factor, polymer cladding / polymer interface')
    plt.title('Weak condition, $F^2$ factor, mirror / polymer interface')
#     np.savetxt('F_factor_Weak_kappa%s.txt' % round(kappa_ext, 3), np.transpose(
#         [(L_ext * 1000), F_factor]), fmt='%1.8e', delimiter='\t')

plt.axvline(x=L_fr, linestyle='dashed', color='b',
            label="$\L_{f_r}$=%s" % round(L_fr, 3))
plt.axvline(x=L_c, linestyle='dashed', color='r',
            label="$\L_{c}$=%s" % round(L_c, 3))
# plt.title('$F^2$ factor, air / polymer interface')
# plt.title('$F^2$ factor, mirror / polymer interface')
plt.xlabel('external cavity length $L$ [mm]')
plt.ylabel('reduction factor $F^2$')
plt.legend(loc=0)

stop = timeit.default_timer()
print('Time: ', stop - start)

plt.show()


# """ 4 Combined plot by function """
# def c_condition:
# def get_color(condition):
#     if condition < 1:
#         return 'red'
#     else:
#         return 'blue'
# F_factor = []
# C_factor = []
# for l_ext in L_ext:
#     tau_EXT = 2 * n_pol * l_ext * 1000 / c
#     X = tau_EXT / tau_L * kappa_ext
#     C = X * np.sqrt(1 + np.square(alpha))
#     C_factor.append(C)
#     print ("C= %s, l_ext= %s" % (C, l_ext))
#     if C > 1 or kappa_ext > 0.5 or R_ext > 0.1:    # strong condition
#         F_ext = R_ext
#         A = tau_EXT / tau_L
#         F = 1 + F_ext * A
#         F_sqr = np.square(F)
#         F_factor.append(F_sqr)
#     elif C < 1 and kappa_ext < 0.5:       # weak condition
#         F = 1 + C
#         F_sqr = np.square(F)
#         F_factor.append(F_sqr)
#     else:
#         F_sqr = 0
#         F_factor.append(F_sqr)
#
# # print len(F_factor)
# plt.figure(1)
# plt.plot(L_ext, F_factor, label="kappa %s, R_ext %s" %
#          (round(kappa_ext, 3), round(R_ext, 3)))
# # plt.axvline(x=L_fr, linestyle='dashed')
# # plt.text(L_fr, 65, r'$\L_{f_2}$=%s' % L_fr, fontsize=15)
# # plt.title('$F^2$ factor, air / polymer interface')
# plt.title('$F^2$ factor, mirror / polymer interface')
# plt.xlabel('external cavity length $L$ [mm]')
# plt.ylabel('reduction factor $F^2$')
# plt.legend(loc=0)
# plt.show()
