#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

E = np.loadtxt("efld_totE.txt")
t = np.loadtxt("efld_time.txt")
plt.plot(t, E)

E_test = np.loadtxt("efld_totE_test.txt")
t_test = np.loadtxt("efld_time_test.txt")
plt.plot(t_test, E_test)


plt.show()
