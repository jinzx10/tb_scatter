import numpy as np
import matplotlib.pyplot as plt

data10 = np.loadtxt("gamma_10.txt")
data4 = np.loadtxt("gamma_4.txt")
data2 = np.loadtxt("gamma_2.txt")
data1 = np.loadtxt("gamma_1.txt")

z10 = data10[:,0]
gamma10 = data10[:,1] / 100
z4 = data4[:,0]
gamma4 = data4[:,1] / 16
z2 = data2[:,0]
gamma2 = data2[:,1] / 4
z1 = data1[:,0]
gamma1 = data1[:,1]

print(np.sum(gamma1))
print(np.sum(gamma2))
print(np.sum(gamma4))
print(np.sum(gamma10))

plt.plot(z1,gamma1)
plt.plot(z2,gamma2)
plt.plot(z4,gamma4)
plt.plot(z10,gamma10)
plt.show()
