#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def lor(x, x0, width):
    return 1/np.pi * (width/2) / ((x-x0)**2 + width**2)

Egrid = np.linspace(0, 10, 1000)
E = np.linspace(0, 10, 20)
broaden_width = 0.5

dos = np.zeros(1000)

for i in E:
    dos = dos + lor(Egrid, i, broaden_width)

plt.plot(Egrid, dos)
plt.show()


