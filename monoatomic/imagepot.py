import numpy as np
import matplotlib.pyplot as plt

z = np.linspace(0, 8, 200)
D = 0.19
C = 2.34
z0 = 0
V = -D / np.sqrt(C*C+(z-z0)**2)

plt.plot(z,V)
plt.show()
