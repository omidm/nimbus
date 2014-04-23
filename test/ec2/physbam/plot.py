#!/usr/bin/env python

import sys
import os
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from decimal import *
getcontext().prec = 2

file_name = sys.argv[1]
f = open(file_name, 'r')
content = f.readlines()

D = []
E = []
T = []
n_plots = 0;

collect = False
for line in content:
  if not collect:
    if line == 'Start\n':
      collect = True;
      d = [];
      e = [];
  else:
    if line == 'End\n':
      collect = False
      D.append(d)
      E.append(e)
      T.append(t)
      n_plots += 1
      continue

    result = re.findall('(\w+):\s*(.*)', line)
    tag = result[0][0]
    data = result[0][1]
    if tag == 'd':
      d = data.split()
    if tag == 'e':
      e = data.split()
    if tag == 't':
      t = data
    
    
    

n_groups = len(D[0])
index = np.arange(n_groups)
y_ticks = 3
font = FontProperties()
font.set_family('serif')

fig, AX = plt.subplots(nrows=n_plots, sharex=True)

print len(AX)
idx = 0;
for idx in np.arange(len(AX)):

  data = np.array(map(float, D[idx]))
  error = np.array(map(float, E[idx]))

  AX[idx].errorbar(index, data, yerr=error, fmt='o')
  AX[idx].set_title(T[idx], fontproperties=font)
  # AX[idx].set_yscale('log')
  # AX[idx].set_ylabel('log time(ms)', fontproperties=font)
  AX[idx].locator_params(axis = 'y', nbins = y_ticks)
  AX[idx].set_ylabel('time(ms)', fontproperties=font)

  idx += 1


plt.xlim((-1, n_groups))
plt.xlabel('Partition', fontproperties=font)
plt.xticks(index, ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII'), fontproperties=font)

# plt.grid()
plt.savefig("test.png")
plt.show()


