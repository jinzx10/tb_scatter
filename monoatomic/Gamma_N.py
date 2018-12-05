#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt("data_1.txt")
data2 = np.loadtxt("data_2.txt")
data4 = np.loadtxt("data_4.txt")
data6 = np.loadtxt("data_6.txt")
data10 = np.loadtxt("data_10.txt")

E = data1[:,0]

Gamma1 = data1[:,1]
dos1 = data1[:,2]

Gamma2 = data2[:,1]
dos2 = data2[:,2]

Gamma4 = data4[:,1]
dos4 = data4[:,2]

Gamma6 = data6[:,1]
dos6 = data6[:,2]

Gamma10 = data10[:,1]
dos10 = data10[:,2]

plt.plot(E, Gamma1)
#plt.plot(E, Gamma2)
plt.plot(E, Gamma4)
#plt.plot(E, Gamma6)
plt.plot(E, Gamma10)


plt.show()
