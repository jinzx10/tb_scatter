#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

t = np.loadtxt("time.txt")
E = np.loadtxt("energy.txt")
#z = np.loadtxt("zcoor.txt")
#q = np.loadtxt("charge.txt")
#vz = np.loadtxt("zvelo.txt")

t_var = np.loadtxt("time_test.txt")
E_var = np.loadtxt("energy_test.txt")

plt.plot(t, E, c='b')
plt.plot(t_var, E_var, c='r')
#plt.plot(t, E, c = 'b', t_var, E_var, c = 'r')

plt.show()
