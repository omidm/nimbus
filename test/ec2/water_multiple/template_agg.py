#!/usr/bin/env python

import decimal
import argparse
import os
import os.path
import re

import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from operator import add
import operator


## Parse the command line arguments ##
parser = argparse.ArgumentParser(description='Plot progress plots.')
parser.add_argument(
    "-s", "--scale",
    dest="scale",
    default=256,
    help="scale, ether 256 or 512")

args = parser.parse_args()

S = int(args.scale)

# 256
idle_time_512 = [10.63, 5.41, 5.03, 3.75, 3.07, 2.6]
idle_time_512 = [100 * t / (t + 25.09) for t in idle_time_512]

idle_time_256 = [1.5 * t for t in idle_time_512]

# Plot the results in stack bar 
N = len(idle_time_256)
ind = np.arange(N)
width = 0.5
# width = 0.25

Legends = []
Legends.append('Scale 512, 64 partitions')
Legends.append('Scale 256, 64 partitions')

n_colors = 3
color_map = cmx.ScalarMappable(
    colors.Normalize(vmin=0, vmax=n_colors+1), cmap='Greens')
Colors = [color_map.to_rgba(i) for i in range(n_colors, -1, -1)]

# Colors = []
# Colors.append('#1CAC78')
# Colors.append('#76FF7A')

Parts = []

# p = plt.bar(map(add, ind, [0.1255] * N), idle_time_512, width, color=Colors[2])
p = plt.bar(ind, idle_time_512, width, color=Colors[2])
for num, rect in zip(idle_time_512, p):
  plt.text(rect.get_x() + rect.get_width()/2.,
           rect.get_y() + rect.get_height()/2.,
           '{:.1f}'.format(num) + '%',
           ha='center', va='center', fontsize=30)

Parts.append(p[0])

# p = plt.bar(map(add, ind, [-0.1255] * N), idle_time_256, width, color=Colors[0])
# for num, rect in zip(idle_time_256, p):
#   plt.text(rect.get_x() + rect.get_width()/2.,
#            rect.get_y() + rect.get_height()/2.,
#            '{:.1f}'.format(num) + '%',
#            ha='center', va='center')
# 
# Parts.append(p[0])



plt.ylabel('Controller Overhead (%)', size=40)

# 256
xticks = ['No Templates', '+ Controller\nTemplate', '+ Version\nMemoization',
       '+ Binding\nMemoization', '+ DM Query\nCache', '+ Worker\nTemplate']

plt.xticks(np.arange(N) + width/2., xticks, fontsize=25)
plt.yticks(fontsize=30)
plt.xlim(0 - width/2., N - width/2.)
# plt.xlim(0 - 2 * width, N - width)

plt.legend(Parts, Legends)
plt.legend(Parts, Legends, loc='upper right', prop={'size':40})

plt.savefig("test.png")
plt.show()


