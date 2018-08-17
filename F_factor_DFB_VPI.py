'''
Created on 15.01.2018

@author: tqian
'''

import scipy
import numpy as np
import matplotlib.pyplot as plt


n_gain = 3.7
n_wg = 3.7
n_air = 1
L_gain = 350
L_ext = np.arange(0, 4, 0.05)  # unit in mm
# L_ext = 3
alpha = 3
c = 3e14    # unit in um
# r_ext = (n_pol - n_air) / (n_pol + n_air)   # 0.183673469388
r_ext = 0.15
R_ext = np.square(r_ext)
print ("R_ext= %s" % R_ext)
R_bragg = 0.3
r_bragg = np.sqrt(R_bragg)
tau_ext = 2 * n_wg * L_ext * 1000 / c
tau_L = (2 * n_gain * L_gain) / c

ratio = R_ext / R_bragg
print ("R_ext/R_bragg= %s" % ratio)

kappa_ext = float(r_ext / r_bragg * (1 - np.square(r_bragg)))
kappa = 1 / tau_L * kappa_ext
print ("kappa_ext= %s" % kappa_ext)
print (kappa)

""" Combined plot """
F_factor = []
for l_ext in L_ext:
    tau_EXT = 2 * n_wg * l_ext * 1000 / c
    X = tau_EXT / tau_L * kappa_ext
    C = X * np.sqrt(1 + np.square(alpha))
    print ("C= %s, l_ext= %s" % (C, l_ext))
    if C > 1 or (kappa_ext > 0.1 or R_ext > 0.1):    # strong condition
        F_ext = R_ext
        A = tau_EXT / tau_L
        F = 1 + F_ext * A
        F_sqr = np.square(F)
        F_factor.append(F_sqr)
    elif C < 1 and kappa_ext < 0.1:       # weak condition
        F = 1 + C
        F_sqr = np.square(F)
        F_factor.append(F_sqr)
    else:
        F_sqr = 0
        F_factor.append(F_sqr)
plt.figure(1)
plt.plot(L_ext, F_factor, label="kappa %s, R_ext %s" %
         (round(kappa_ext, 3), round(R_ext, 3)))
# plt.title('$F^2$ factor, air / polymer interface')
plt.title('$F^2$ factor, mirror / polymer interface')
plt.xlabel('external cavity length $L$ [mm]')
plt.ylabel('reduction factor $F^2$')
plt.legend(loc=0)
plt.show()
