#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from mpl_toolkits.mplot3d import Axes3D

num_frames = 1000;

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

x_lat = np.loadtxt("efld_lat_coor_x.txt")
y_lat = np.loadtxt("efld_lat_coor_y.txt")
z_lat = np.loadtxt("efld_lat_coor_z.txt")

x_mol = np.loadtxt("efld_mol_coor_x.txt")
y_mol = np.loadtxt("efld_mol_coor_y.txt")
z_mol = np.loadtxt("efld_mol_coor_z.txt")

x = np.concatenate([x_lat, x_mol], axis=1)
y = np.concatenate([y_lat, y_mol], axis=1)
z = np.concatenate([z_lat, z_mol], axis=1)

n_raw = len(x)
spacing = int(n_raw/num_frames)
print(spacing)
x = x[::spacing]
y = y[::spacing]
z = z[::spacing]

scat = ax.scatter(x[0], y[0], z[0])
#plt.show()

def update(i):
    scat._offsets3d = (x[i], y[i], z[i])

ani = ani.FuncAnimation(fig, update, np.arange(1,num_frames), interval=25, blit=False)
ani.save('animation.mp4')

