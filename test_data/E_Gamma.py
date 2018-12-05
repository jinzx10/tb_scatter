#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

E_Gamma_nk16 = np.loadtxt("E_Gamma_z4_nk16.txt")
E_Gamma_nk36 = np.loadtxt("E_Gamma_z4_nk16.txt")

E_16 = E_Gamma_nk16[:, 0]
E_36 = E_Gamma_nk36[:, 0]

Gamma_on_top_nk16 = E_Gamma_nk16[:, 1]
Gamma_bridge_nk16 = E_Gamma_nk16[:, 1]
Gamma_hollow_nk16 = E_Gamma_nk16[:, 2]

Gamma_on_top_nk36 = E_Gamma_nk36[:, 1]
Gamma_bridge_nk36 = E_Gamma_nk36[:, 1]
Gamma_hollow_nk36 = E_Gamma_nk36[:, 2]

plt.plot(E_16, Gamma_on_top_nk16)
plt.plot(E_16, Gamma_bridge_nk16)
plt.plot(E_16, Gamma_hollow_nk16)

#plt.plot(E_36, Gamma_on_top_nk36, marker='o', linestyle='None')
#plt.plot(E_36, Gamma_bridge_nk36, marker='o', linestyle='None')
#plt.plot(E_36, Gamma_hollow_nk36, marker='o', linestyle='None')

plt.show()
