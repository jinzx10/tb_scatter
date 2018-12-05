#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("E_dep.txt")
E = data[:, 0]
Gamma = data[:, 1]
dos = data[:, 2]

plt.plot(E, Gamma)
#plt.plot(E, dos)

plt.show()
