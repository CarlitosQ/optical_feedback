'''
Created on 02.07.2018

@author: tqian
'''

import scipy
import numpy as np
from numpy import pi
import cmath
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

c = 3e8
freq = np.linspace(0, 100, 10000)
freq = freq * 1e9
hv = 1.5
a = 5e-16 * 1e-4     # 5e-16 cm^2
tau = 3e-9
tau_p = 2e-12
yeta_i = 0.867
yeta_d = 0.8
v_g = 3e10 / 4 * 1e-2
V_p = 5 * 0.25 * 200e-18
alpha_i = 5e2
alpha_m = 60e2
q = 1.6e-19
I = 100e-3
I_th = 20e-3
g_th = alpha_i + alpha_m
g_th = 1 / (v_g * tau_p)

N_p0 = yeta_i * (I - I_th) / (q * v_g * g_th * V_p)
omega = 2 * pi * freq
omega_r = np.sqrt(v_g * a * N_p0 / tau_p)
modulation = (yeta_d * hv / q) / (1 - (omega / omega_r)**2 + 1j * (omega / omega_r)
                                  * (omega_r * tau_p + 1 / omega_r * tau))

plt.figure(1)
plt.plot(freq, np.log(modulation))
plt.show()