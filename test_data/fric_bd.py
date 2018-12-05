#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

fric_bd16 = np.loadtxt("fric_bd16.txt")
width16 = fric_bd16[:, 0]
fric16 = fric_bd16[:, 1]
plt.plot(width16, fric16)

fric_bd36 = np.loadtxt("fric_bd36.txt")
width36 = fric_bd36[:, 0]
fric36 = fric_bd36[:, 1]
plt.plot(width36, fric36)

#fric_bd100 = np.loadtxt("fric_bd100.txt")
#width100 = fric_bd100[:, 0]
#fric100 = fric_bd100[:, 1]
#plt.plot(width100, fric100)
#
#fric_bd400 = np.loadtxt("fric_bd400.txt")
#width400 = fric_bd400[:, 0]
#fric400 = fric_bd400[:, 1]
#plt.plot(width400, fric400)

#fric_bd900 = np.loadtxt("fric_bd900.txt")
#width900 = fric_bd900[:, 0]
#fric900 = fric_bd900[:, 1]
#plt.plot(width900, fric900)

plt.show()

