#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt("cpl_pot.txt")
x = data[0,:]
y = data[1,:]
V = data[2,:]


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
graph = ax.scatter(x, y, V)
plt.show()
