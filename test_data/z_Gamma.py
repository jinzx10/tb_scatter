#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

z_Gamma_nk16 = np.loadtxt("z_Gamma_nk16.txt")
z_Gamma_nk36 = np.loadtxt("z_Gamma_nk36.txt")

z16 = z_Gamma_nk16[:,0]
z36 = z_Gamma_nk36[:,0]

Gamma_ontop_nk16 = z_Gamma_nk16[:,1]
Gamma_bridge_nk16 = z_Gamma_nk16[:,2]
Gamma_hollow_nk16 = z_Gamma_nk16[:,3]

Gamma_ontop_nk36 = z_Gamma_nk36[:,1]
Gamma_bridge_nk36 = z_Gamma_nk36[:,2]
Gamma_hollow_nk36 = z_Gamma_nk36[:,3]

plt.plot(z16, Gamma_ontop_nk16)
plt.plot(z16, Gamma_bridge_nk16)
plt.plot(z16, Gamma_hollow_nk16)

#plt.plot(z36, Gamma_ontop_nk36, marker='o', linestyle='None')
#plt.plot(z36, Gamma_bridge_nk36, marker='o', linestyle='None')
#plt.plot(z36, Gamma_hollow_nk36, marker='o', linestyle='None')


plt.show()
