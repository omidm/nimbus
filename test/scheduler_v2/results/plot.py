#!/usr/bin/env python

import sys
import os
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from decimal import *
getcontext().prec = 2

P = 4

X = (128, 256, 512, 1024)

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

I = (88 - 142, 150 - 203, 211 - 262, 271 - 321, 330 - 382, 393 - 449)

A = (1 * 128 / 10, 2.4 * 256 / 10, 5.1 * 512 / 10, 15.68 * 1024 / 10)

J = (7.57, 22.73, 79.80, 287.78)

a = plt.errorbar(X, A, fmt='s-', color='r', ms=8)
j = plt.errorbar(X, J, fmt='^-', color='b', ms=8)
v = plt.errorbar(X, V, fmt='o-', color='g', ms=8)
hline = plt.axhline(-np.mean(I), color='#fb8072', ls='--', linewidth=4)

plt.legend((a, j, v, hline),
    ('Assigning First Stage', 'Spawning Jobs', 'Versioning Job Graph', 'Iteration Duration (128 cube)'),
    loc='upper left')
plt.yscale('log')
plt.xscale('log')
plt.xlim(100, 1025)

plt.ylabel('Time (seconds)')
plt.xlabel('Number of Application Partitions')
plt.title('Analyzing Scalability of Nimbus in Water Simulation')

# plt.grid()
plt.savefig("test.png")
plt.show()


