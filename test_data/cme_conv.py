#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

t20 = np.loadtxt("time_dt20.txt")
E20 = np.loadtxt("energy_dt20.txt")
#z = np.loadtxt("zcoor.txt")
#q = np.loadtxt("charge.txt")
#vz = np.loadtxt("zvelo.txt")

t10 = np.loadtxt("time_dt10.txt")
E10 = np.loadtxt("energy_dt10.txt")

plt.plot(t10, E10, c='b')
plt.plot(t20, E20, c='r')
#plt.plot(t, E, c = 'b', t_var, E_var, c = 'r')

plt.show()
