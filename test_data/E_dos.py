#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

E_dos = np.loadtxt("E_dos_z4_nk16.txt")
E = E_dos[:, 0]
dos = E_dos[:, 1]


plt.plot(E, dos)

plt.show()
