import numpy as np
import matplotlib.pyplot as plt
import sys

crystal = float(sys.argv[2])

data = np.loadtxt('DirectAbs.dat')
w = data[:, 0]
a = data[:, 1]
a /= np.max(a)
plt.plot(w + 19000 - crystal, a, label='Ian')

try:
    data = np.loadtxt(sys.argv[1])
    w = data[:, 0]
    a = data[:, 1]
    plt.plot(w, a, label='April')
except:
    pass

plt.legend()
plt.title('Absorption Spectra for TP-PDI')
plt.xlabel('$\omega\;(cm^{-1})$')
plt.ylabel('$A(\omega)$')
plt.show()
