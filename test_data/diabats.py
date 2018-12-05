#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

diabats = np.loadtxt("diabats.txt")
diabats_var = np.loadtxt("diabats_var.txt")
adiabats = np.loadtxt("adiabats.txt")
bd_diabats = np.loadtxt("broadened_diabats.txt")

z = diabats[:, 0]
z_var = diabats_var[:, 0]

# diabats
neutral_on_top = diabats[:, 1]
neutral_bridge = diabats[:, 2]
neutral_hollow = diabats[:, 3]

ionic_on_top = diabats[:, 4]
ionic_bridge = diabats[:, 5]
ionic_hollow = diabats[:, 6]

# diabats_var
neutral_on_top_var = diabats_var[:, 1]
neutral_bridge_var = diabats_var[:, 2]
neutral_hollow_var = diabats_var[:, 3]

ionic_on_top_var = diabats_var[:, 4]
ionic_bridge_var = diabats_var[:, 5]
ionic_hollow_var = diabats_var[:, 6]

# broadened diabats
bd_neutral_on_top = bd_diabats[:, 1]
bd_neutral_bridge = bd_diabats[:, 2]
bd_neutral_hollow = bd_diabats[:, 3]

bd_ionic_on_top = bd_diabats[:, 4]
bd_ionic_bridge = bd_diabats[:, 5]
bd_ionic_hollow = bd_diabats[:, 6]

# adiabats
adiabats_on_top = adiabats[:, 1]
adiabats_bridge = adiabats[:, 2]
adiabats_hollow = adiabats[:, 3]


# plot
#plt.plot(z, neutral_on_top, color = (0.8, 0.8, 1), marker = 'None', linestyle='solid')
#plt.plot(z, ionic_on_top, color = (1, 0.8, 0.8), marker = 'None', linestyle='solid')
#plt.plot(z, adiabats_on_top, color = (0.8, 1, 0.8), marker = 'None', linestyle='solid')
#plt.plot(z, bd_neutral_on_top, color = (0.8, 0.8, 1), marker = 'o', linestyle='None', markerfacecolor='None')
#plt.plot(z, bd_ionic_on_top, color = (1, 0.8, 0.8), marker = 'o', linestyle='None', markerfacecolor='None')
#plt.plot(z_var, neutral_on_top_var, color=(0.8, 0.8, 1), marker = 'o', linestyle='None', markerfacecolor='None')
#plt.plot(z_var, ionic_on_top_var, color=(1, 0.8, 0.8), marker = 'o', linestyle='None', markerfacecolor='None')

#plt.plot(z, neutral_bridge, color = (0.4, 0.4, 1), marker = 'None', linestyle='solid')
#plt.plot(z, ionic_bridge, color = (1, 0.4, 0.4), marker = 'None', linestyle='solid')
#plt.plot(z, adiabats_bridge, color = (0.4, 1, 0.4), marker = 'None', linestyle='solid')
#plt.plot(z, bd_neutral_bridge, color = (0.4, 0.4, 1), marker = 'o', linestyle='None', markerfacecolor='None')
#plt.plot(z, bd_ionic_bridge, color = (1, 0.4, 0.4), marker = 'o', linestyle='None', markerfacecolor='None')
#plt.plot(z_var, neutral_bridge_var, color=(0.4, 0.4, 1), marker = 'o', linestyle='None', markerfacecolor='None')
#plt.plot(z_var, ionic_bridge_var, color=(1, 0.4, 0.4), marker = 'o', linestyle='None', markerfacecolor='None')

plt.plot(z, neutral_hollow, color = (0, 0, 1), marker = 'None', linestyle='solid')
plt.plot(z, ionic_hollow, color = (1, 0, 0), marker = 'None', linestyle='solid')
plt.plot(z, adiabats_hollow, color = (0, 1, 0), marker = 'None', linestyle='solid')
plt.plot(z, bd_neutral_hollow, color = (0, 0, 1), marker = 'o', linestyle='None', markerfacecolor='None')
plt.plot(z, bd_ionic_hollow, color = (1, 0, 0), marker = 'o', linestyle='None', markerfacecolor='None')
#plt.plot(z_var, neutral_hollow_var, color=(0, 0, 1), marker = 'o', linestyle='None', markerfacecolor='None')
#plt.plot(z_var, ionic_hollow_var, color=(1, 0, 0), marker = 'o', linestyle='None', markerfacecolor='None')



plt.show()

