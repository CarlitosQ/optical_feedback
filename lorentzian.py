# Lorentzian lineshape

import scipy
import numpy as np
import matplotlib.pyplot as plt


def lorentzian_shape(tau, x, x_0):
    return (1/np.pi)*((1/2*tau)/((x-x_0)**2+(1/2*tau)**2))


x = np.linspace(-10, 10, 1000)
tau = 5
lorentz = lorentzian_shape(tau, x, 0)
# x1 = np.linspace(-6, , 1000)
# x2 = np.linspace(-3, 6, 1000)
plt.figure(1)
plt.plot(x-11, lorentz, color='black')
plt.plot(x, lorentz, color='black')
plt.axis('off')

plt.figure(2)
lorentz = lorentzian_shape(6.5, x, 0)
plt.plot(x, lorentz, color='black')
plt.axis('off')
plt.show()
