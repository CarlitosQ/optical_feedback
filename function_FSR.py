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

L_grat = 699.4075 * 1e-6
L_grat_eff = L_grat / 2
L_poly = (180 + 190 + 155) * 1e-6 + L_grat_eff
L_wg = (180 + 190 + 155) * 1e-6

r1 = np.sqrt(0.9)
r_1 = 0.99
r_ext = 0.5
n_gain = 3.2
n_air = 1.
n_poly = 1.46

def FreeSpectralRange(L_gain, L_ext):
    FreeSpectralRange = (c / (2 * (n_gain * L_gain*1e-6 + n_poly * (L_poly + L_ext*1e-6))))/1e9
    return FreeSpectralRange

FSR = FreeSpectralRange(300, 0)
FSR = FreeSpectralRange(500, 0)
print ("FSR = %s GHz" % float('%.4g' % FSR))