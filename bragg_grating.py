'''
Created on Mar 9, 2018

@author: tqian
'''

import scipy
import numpy as np
import matplotlib.pyplot as plt

# kappa_dc=n*omega*ips*
L = np.arange(0, 300, 0.05)
kappa_dc = -31.91 - 75.33j
kappa_ac = kappa_dc / 2
delta = kappa_dc
alpha = np.sqrt(np.square(abs(kappa_ac)) - np.square(delta))
r = (-kappa_ac * np.sinh(alpha * L)) / \
    (delta * np.sinh(alpha * L) - 1j * alpha * np.cosh(alpha * L))
plt.figure(1)
plt.plot(r)
plt.show()
