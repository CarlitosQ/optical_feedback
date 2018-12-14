import scipy
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal


def F_factor(L_grat, L_wg_before_grat, L_ext, r_ext):
    c = 3e8
    n_gain = 3.2
    n_poly = 1.46
    R_bragg = 0.3
    r_bragg = np.sqrt(R_bragg)
    L_gain = 300e-6
    L_grat = L_grat
    L_grat_eff = L_grat / 2
    L_poly = (L_wg_before_grat + L_grat_eff) * 1e-6    # 170 phase section
    alpha = -3

    tau_ext = 2 * n_poly * L_ext / c
    tau_L = 2 * (n_gain * L_gain + n_poly * L_poly) / c

    # 1. practical equation
    yeta_couple = 1
    R_ext = r_ext**2
    F_ext = R_ext * yeta_couple
    A = tau_ext / tau_L
    F = 1 + F_ext * A

    # 2. theretical equation
    # kappa_ext = float(r_ext / r_bragg * (1 - np.square(r_bragg)))
    # X = tau_ext / tau_L * kappa_ext
    # C = X * np.sqrt(1 + alpha**2)
    # F = 1 + C

    F_sqr = np.square(F)
    return F_sqr


def C_factor(L_grat, L_wg_before_grat, L_ext, r_ext):
    c = 3e8
    n_gain = 3.2
    n_poly = 1.46
    R_bragg = 0.3
    r_bragg = np.sqrt(R_bragg)
    L_gain = 300e-6
    L_grat = L_grat
    L_grat_eff = L_grat / 2
    L_poly = (L_wg_before_grat + L_grat_eff) * 1e-6    # 170 phase section
    alpha = -3

    tau_ext = 2 * n_poly * L_ext / c
    tau_L = 2 * (n_gain * L_gain + n_poly * L_poly) / c

    # 1. practical equation
    # yeta_couple = 1
    # R_ext = r_ext**2
    # F_ext = R_ext * yeta_couple
    # A = tau_ext / tau_L
    # F = 1 + F_ext * A

    # 2. theretical equation
    kappa_ext = float(r_ext / r_bragg * (1 - np.square(r_bragg)))
    X = tau_ext / tau_L * kappa_ext
    C = X * np.sqrt(1 + alpha**2)

    return C


def normal_cavity(L_grat, L_wg_before_grat):
    c = 3e8
    n_gain = 3.2
    n_poly = 1.46
    R_bragg = 0.3
    r_bragg = np.sqrt(R_bragg)
    L_gain = 300e-6
    L_grat = L_grat
    L_grat_eff = L_grat / 2
    L_poly = (L_wg_before_grat + L_grat_eff)    # 170 phase section

    normal_cavity = c / (2 * (n_gain * L_gain + n_poly * L_poly))
    return normal_cavity


def external_cavity(L_grat, L_wg_before_grat, L_ext, r_ext):
    c = 3e8
    n_gain = 3.2
    n_poly = 1.46
    R_bragg = 0.3
    r_bragg = np.sqrt(R_bragg)
    L_gain = 300e-6
    L_grat = L_grat
    L_grat_eff = L_grat / 2
    L_poly = (L_wg_before_grat + L_grat_eff)    # 170 phase section

    external_cavity = c / (2 * n_poly * (L_ext + L_grat_eff))
    return external_cavity


def FSR(L_grat, L_wg_before_grat, L_ext, r_ext):
    c = 3e8
    n_gain = 3.2
    n_poly = 1.46
    R_bragg = 0.3
    r_bragg = np.sqrt(R_bragg)
    L_gain = 300e-6
    L_grat = L_grat
    L_grat_eff = L_grat / 2
    L_poly = (L_wg_before_grat + L_grat_eff)   # 170 phase section

    FSR = c / (2 * (n_gain * L_gain + n_poly *
                    (L_poly + (L_ext + L_grat_eff))))
    return FSR


# Set 1
# n_poly = 1.46
# n_air = 1.
# r_ext = (n_poly - n_air) / (n_poly + n_air)   # 0.18699186991869918
# # L_grat = 699.84 * 1e-6
# L_grat = 700.644*1e-6
# L_grat_eff = L_grat / 2
# # L_wg_before_grat = 525e-6
# L_wg_before_grat = 509e-6
# # L_ext = (267 + 518.362787842 + 688) * 1e-6 + L_grat_eff
# L_ext = 6359.758538 * 1e-6
# # print(F_factor(L_grat, L_wg_before_grat, L_ext, r_ext))
# # print(C_factor(L_grat, L_wg_before_grat, L_ext, r_ext))

# Set 2
n_poly = 1.46
n_air = 1.
r_ext = (n_poly - n_air) / (n_poly + n_air)   # 0.18699186991869918
L_grat = 699.84 * 1e-6
L_grat_eff = L_grat / 2
L_wg_before_grat = 525e-6
L_ext = (267 + 518.362787842 + 688) * 1e-6
# print(F_factor(L_grat, L_wg_before_grat, L_ext, r_ext))
print(C_factor(L_grat, L_wg_before_grat, L_ext, r_ext))

print("mode spacing for ADVA = %s GHz" %
      (normal_cavity(L_grat, L_wg_before_grat) / 1e9))
print("external cavity spacing = %s GHz" %
      (external_cavity(L_grat, L_wg_before_grat, L_ext, r_ext) / 1e9))
print("FSR = %s GHz" % (FSR(L_grat, L_wg_before_grat, L_ext, r_ext)/1e9))
