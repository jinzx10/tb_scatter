#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

coor = np.loadtxt("k.txt")
band = np.loadtxt("bandE.txt")
kx = coor[0,:]
ky = coor[1,:]
#z = coor[2,:]
#b2 = band[1,:]


graph = ax.scatter(kx,ky,band)
#graph = ax.scatter(kx,ky,b2)
#graph = ax.scatter(1,2,3)


plt.show()
