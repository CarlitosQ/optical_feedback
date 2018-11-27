import scipy
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal

c = 3e8
n_gain = 3.2
n_poly = 1.46
n_air = 1.
L_grat = 699.84 * 1e-6
L_grat_eff = L_grat / 2
L_poly = (525 + L_grat_eff) * 1e-6    # 170 phase section
L_gain = 300e-6
alpha = -3

x_r = np.sqrt(L_gain*np.real((1+alpha*1j)/(L_gain+L_grat)))
print(x_r)
