#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("z_dep.txt")
z = data[:,0]

neutral_ontop = data[:,1]
ionic_ontop = data[:,2]

neutral_hollow = data[:,3]
ionic_hollow = data[:,4]

neutral_bridge = data[:,5]
ionic_bridge = data[:,6]

neutral_ehollow = data[:,7]
ionic_ehollow = data[:,8]

#plt.plot(z,neutral_ontop)
#plt.plot(z,ionic_ontop)
plt.plot(z,neutral_hollow)
plt.plot(z,ionic_hollow)
plt.plot(z,neutral_bridge)
plt.plot(z,ionic_bridge)
plt.plot(z,neutral_ehollow)
plt.plot(z,ionic_ehollow)

plt.show()
