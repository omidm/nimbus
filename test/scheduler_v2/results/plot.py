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


# collected data
I = (88 - 142, 150 - 203, 211 - 262, 271 - 321, 330 - 382, 393 - 449)

VV = []
VV.append((0.68, 0.61, 0.66, 0.67, 0.57))
VV.append((2.02, 1.76, 2.09, 2.01, 1.96))
VV.append((7.17, 7.16, 7.25, 8.20, 8.35))
VV.append((27.77, 27.35, 27.22, 27.71, 27.34))
V = []
E = []
for i in range(0, len(VV)):
  V.append(np.mean(VV[i]))
  E.append(np.std(VV[i]))

AO = (1 * 128 / 10, 2.4 * 256 / 10, 5.1 * 512 / 10, 15.68 * 1024 / 10)

AN = (0.05 * 128 / 10, 0.043 * 256 / 10, 0.041 * 512 / 10, 0.044 * 1024 / 10)

J = (7.57, 22.73, 79.80, 287.78)


# prepare data for plotting
X = (128, 256, 512, 1024)

Y = []
Y.append(AO)
Y.append(AN)
Y.append(J)
Y.append(V)

Colors = []
Colors.append('r')
Colors.append('r')
Colors.append('b')
Colors.append('g')

Formats = []
Formats.append('s-')
Formats.append('s--')
Formats.append('^-')
Formats.append('o-')

Legends = []
Legends.append('Assigning First Stage (Old)')
Legends.append('Assigning First Stage (New)')
Legends.append('Spawning Jobs')
Legends.append('Versioning Job Graph')

P = []

for i in range(0, len(Y)):
  P.append(plt.errorbar(X, Y[i], fmt=Formats[i], color=Colors[i], ms=8))
  for j in range(0, N):
    tag = '%2.2f' % Y[i][j]
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


