#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

kin_loss = np.loadtxt("kin_loss.txt")
out_angle = np.loadtxt("out_angle.txt")

kin_loss_percent = kin_loss[:,2]

#plt.hist(kin_loss_percent)
plt.hist(out_angle, bins=20)
plt.show()
