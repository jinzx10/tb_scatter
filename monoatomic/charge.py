#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

t = np.loadtxt("time.txt")
q = np.loadtxt("charge.txt")

#plt.plot(t,q)
plt.plot(q)
plt.show()
