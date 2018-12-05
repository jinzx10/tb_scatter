import numpy as np
import matplotlib.pyplot as plt

r0 = 5.5
l = 0.4
C = 0.03
r = np.linspace(3, 8, 100)

t = C / (1 + np.exp( (r-r0)/l ))

plt.plot(r,t)
plt.show()
