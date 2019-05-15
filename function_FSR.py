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
L_wg = (180 + 190 + 155) * 1e-6
L_poly = L_wg + L_grat_eff


r1 = np.sqrt(0.9)
r_1 = 0.99
r_ext = 0.5
n_gain = 3.2
n_air = 1.
n_poly = 1.46


def FreeSpectralRange(L_gain, L_grat, L_ext):
    L_wg = (180 + 190 + 155) * 1e-6
    L_poly = L_wg + L_grat
    FSR = (c / (2 * (n_gain * L_gain * 1e-6 + n_poly * (L_poly + L_ext)))) / 1e9
    return FSR


L_grat_700 = 700.644*1e-6
L_grat_300 = 299.7*1e-6
L_ext_700 = (266.45 + 518.362787842 + 688) * 1e-6  # Lgrating 700
L_ext_300 = (264.49 + 402.17807 + 518.362787842 + 688)*1e-6  # Lgrating 300
FSR_700 = FreeSpectralRange(300, L_grat_700, L_ext_700)
FSR_300 = FreeSpectralRange(300, L_grat_300, L_ext_300)
# FSR = FreeSpectralRange(500, 0)
print("FSR = %s GHz" % float('%.4g' % FSR_700))
print("FSR = %s GHz" % float('%.4g' % FSR_300))
