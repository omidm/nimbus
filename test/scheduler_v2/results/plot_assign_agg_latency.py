#!/usr/bin/env python

import sys
import os
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from decimal import *
getcontext().prec = 2

N = 4


I = (88 - 142, 150 - 203, 211 - 262, 271 - 321, 330 - 382, 393 - 449)

C = (1.70, 4.03, 10.05, 36.22) 

L = (0.007, 0.012, 0.020, 0.081)

# prepare data for plotting
X = (128, 256, 512, 1024)

Y = []
Y.append(C)
Y.append(L)

Colors = []
Colors.append('r')
Colors.append('g')

Formats = []
Formats.append('s-')
Formats.append('s--')
Formats.append('^-')
Formats.append('o-')

Legends = []
Legends.append('Latency of First Assignment (Complete Versioning)')
Legends.append('Latency of First Assignment (Lazy Versioning)')

P = []

for i in range(0, len(Y)):
  P.append(plt.errorbar(X, Y[i], fmt=Formats[i], color=Colors[i], ms=8))
  for j in range(0, N):
    tag = '%2.3f' % Y[i][j]
    plt.annotate(tag, xy=(X[j], Y[i][j]))


Legends.append('Iteration Duration (128 cube)')
P.append(plt.axhline(-np.mean(I), color='#fb8072', ls='--', linewidth=4))
tag = '%2.2f' % (-np.mean(I))
plt.annotate(tag, xy=(X[0], (-np.mean(I))))

plt.legend(P, Legends, loc='upper left')

plt.yscale('log')
plt.xscale('log')
plt.xlim(100, 1024 * 1.5)

plt.ylabel('Time (seconds)')
plt.xlabel('Number of Application Partitions')
plt.title('Analyzing Scalability of Nimbus in Water Simulation')



# plt.grid()
plt.savefig("test.png")
plt.show()


