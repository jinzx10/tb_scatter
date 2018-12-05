#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

x = np.loadtxt("xcoor.txt")
y = np.loadtxt("ycoor.txt")
z = np.loadtxt("zcoor.txt")
n = len(x)

scat = ax.scatter(x[0], y[0], z[0])
#plt.show()

def update(i):
    scat._offsets3d = (x[i], y[i], z[i])

ani = ani.FuncAnimation(fig, update, np.arange(1,n), interval=25, blit=False)
ani.save('animation.mp4')

