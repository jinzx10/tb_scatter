#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

coor = np.loadtxt("sc_coor.txt")
x = coor[0,:]
y = coor[1,:]
z = coor[2,:]

graph = ax.scatter(x,y,z)
plt.show();
