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
YV = []
YV.append((0.78, 0.71, 0.76, 0.77, 0.67))

Y = (0.68, 2.02, 4.01, 6.0)
# E = Y / 5

plt.errorbar(X, Y, yerr=Y, fmt='o-')


# n_groups = len(D[0])
# index = np.arange(n_groups)
# y_ticks = 5
# font = FontProperties()
# font.set_family('serif')
# 
# fig, AX = plt.subplots(nrows=n_plots, sharex=True)
# 
# print len(AX)
# idx = 0;
# for idx in np.arange(len(AX)):
# 
#   data = np.array(map(float, D[idx]))
#   error = np.array(map(float, E[idx]))
# 
#   AX[idx].errorbar(index, data, yerr=error, fmt='o')
#   AX[idx].set_title(T[idx], fontproperties=font)
#   # AX[idx].set_yscale('log')
#   # AX[idx].set_ylabel('log time(ms)', fontproperties=font)
#   AX[idx].locator_params(axis = 'y', nbins = y_ticks)
#   AX[idx].set_ylabel('time(ms)', fontproperties=font)
#   AX[idx].set_ylim(0, max(data + error) * 1.1)
# 
#   idx += 1
# 
# 
# plt.xlim((-1, n_groups))
# plt.xlabel('Partition', fontproperties=font)
# plt.xticks(index, ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII'), fontproperties=font)

# plt.grid()
plt.savefig("test.png")
plt.show()


