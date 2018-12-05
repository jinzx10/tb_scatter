#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("z_gamma.txt")
z = data[:,0]

gamma_ontop = data[:,1]
gamma_ehollow = data[:,2]
gamma_bridge = data[:,3]
gamma_hollow = data[:,4]

plt.plot(z, gamma_ontop)
plt.plot(z, gamma_ehollow)
plt.plot(z, gamma_bridge)
plt.plot(z, gamma_hollow)

plt.show()
