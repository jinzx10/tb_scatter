#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

k = np.loadtxt("k.txt")
kx = k[0,:]
ky = k[1,:]

plt.scatter(kx,ky)
plt.show()
