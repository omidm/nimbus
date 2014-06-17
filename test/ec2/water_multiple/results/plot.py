#!/usr/bin/env python

import decimal
import argparse
import os
import os.path
import re

import numpy as np
import matplotlib.pyplot as plt
from operator import add

## Parse the command line arguments ##
parser = argparse.ArgumentParser(description='Process log files.')
parser.add_argument(
    "-i", "--input",
    dest="input_base_name",
    default="data",
    help="base input file name that follows with -<worker number>")
parser.add_argument(
    "-d", "--directory",
    dest="dir",
    default=".",
    help="directory to find the input files")
parser.add_argument(
    "-wn", "--workernum",
    dest="worker_num",
    default=1,
    type=int,
    help="number of workers to process")

args = parser.parse_args()

N = args.worker_num
P = 5

# Fetch data from data files
compute_bare = []
compute_trans = []
copy_bare = []
copy_trans = []
idle = []
total = []

for i in range(1, N + 1):
  print 'Processing the data for worker ' + str(i) + '...'

  data_file =  args.input_base_name + '-' + str(i)

  print 'Opening file ' + data_file
  data = open(data_file, 'r')
  
  cmb_found = False
  cmt_found = False
  cob_found = False
  cot_found = False
  i_found   = False
  t_found   = False

  for num, line in enumerate(data):
      x =  re.findall('(\d+\.\d+$|\d+e-\d+$|\d+$)', line)
      if len(x) > 0:
          xnum = decimal.Decimal(x[0])
          if "Application" in line:
              total.append(xnum)
              t_found = True
          if "Decomposed Compute Bare" in line:
              compute_bare.append(xnum)
              cmb_found = True
          if "Decomposed Compute Translator" in line:
              compute_trans.append(xnum)
              cmt_found = True
          if "Decomposed Copy Bare" in line:
              copy_bare.append(xnum)
              cob_found = True
          if "Decomposed Copy Translator" in line:
              copy_trans.append(xnum)
              cot_found = True
          if "Decomposed Idle" in line:
              idle.append(xnum)
              i_found = True

  if (not cmb_found or \
      not cmt_found or \
      not cob_found or \
      not cot_found or \
      not i_found or \
      not t_found):
    print 'File ' + data_file + ' does not match the format.'
    exit(0)

# Plot the results in stack bar 
ind = np.arange(N)
width = 0.35

Data = []
Data.append(compute_bare)
Data.append(compute_trans)
Data.append(copy_bare)
Data.append(copy_trans)
Data.append(idle)

Legends = []
Legends.append('Compute Job')
Legends.append('Translator for Computation')
Legends.append('Copy Job')
Legends.append('Translator for Copy')
Legends.append('Idle')

Colors = []
Colors.append('r')
Colors.append('g')
Colors.append('b')
Colors.append('y')
Colors.append('w')

bottom = []
bottom.append([0] * N)

Parts = []

for i in range(0, P):
  p = plt.bar(ind, Data[i], width, color=Colors[i], bottom=bottom[i])
  bottom.append(map(add, bottom[i], Data[i]))
  Parts.append(p[0])

plt.ylabel('Time (seconds)')
plt.title('Nimbus')
xticks = []
for i in range(1, N + 1):
  xticks.append('W' + str(i))
plt.xticks(ind+width/2., xticks )
plt.yticks(np.arange(0,81,10))
# plt.legend( (S[0][0], S[1][0], S[2][0], S[3][0], S[4][0]), ('Translator Compute', 'Compute', 'Translator Copy', 'Copy', 'Idle') )
plt.legend(Parts, Legends)

plt.show()


