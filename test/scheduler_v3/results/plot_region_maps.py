#!/usr/bin/env python

import sys
import os
import re
import argparse
import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations


def draw_geometric_region(figure, x, y, z, dx, dy, dz, color):
  X = [x, x + dx]
  Y = [y, y + dy]
  Z = [z, z + dz]

  for s, e in combinations(np.array(list(product(X,Y,Z))), 2):
    if np.sum(s == e) == 2:
      ax.plot3D(*zip(s,e), color=color)






## Parse the command line arguments ##
parser = argparse.ArgumentParser(description='Process log files.')
parser.add_argument(
    "-i", "--input",
    dest="ifname",
    default="load_balancer_log",
    help="input file name")
parser.add_argument(
    "-d", "--directory",
    dest="idir",
    default="../",
    help="directory to find the input file")
parser.add_argument(
    "-od", "--outdirectory",
    dest="odir",
    default=".",
    help="directory to dump the output file")
parser.add_argument(
    "-o", "--output",
    dest="ofname",
    default="test",
    help="output file name")

args = parser.parse_args()


Colors = []
Colors.append('k')
Colors.append('r')
Colors.append('b')
Colors.append('g')
Colors.append('y')

file_name = args.idir + '/' + args.ifname
print 'Opening the file ' + file_name
log = open(file_name, 'r')

fig = []
fig_num = 0;
worker_id = 0;

for num, line in enumerate(log):
  if "Region Map Begin" in line:
    fig_num += 1
    fig = plt.figure()
    ax = fig.gca(projection='3d')

  if "worker_id:" in line:
    x =  re.findall('(\d+)', line)
    assert(len(x) == 1)
    worker_id = int(x[0])

  if "bbox:" in line:
    x =  re.findall('(\d+\.\d+|\d+)', line)
    assert(len(x) == 6)
    draw_geometric_region(fig, int(x[0]), int(x[1]), int(x[2]), int(x[3]), int(x[4]), int(x[5]), Colors[worker_id]);

  if "Region Map End" in line:
    plt.savefig(args.ofname + "_" + str(fig_num) + ".png")
    plt.show()
    time.sleep(0.1)
    plt.close()

sys.exit(0)


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")

draw_geometric_region(fig, 0, 0, 0, 1, 2, 4, "b");
draw_geometric_region(fig, -1, -1, -1, 10, 10, 10, "w");

plt.savefig(args.ofname)
plt.show()

sys.exit(0)

