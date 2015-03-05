#!/usr/bin/env python

import decimal
import argparse
import os
import os.path
import re

import numpy as np
import matplotlib.pyplot as plt
from operator import add
import operator

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
    required=True,
    type=int,
    help="number of workers to process")
parser.add_argument(
    "-cn", "--corenum",
    dest="core_num",
    default=1,
    required=True,
    type=int,
    help="number of cores per worker")
parser.add_argument(
    "-in", "--iternum",
    dest="iter_num",
    default=1,
    required=True,
    type=int,
    help="number of iterations")

args = parser.parse_args()

WN = args.worker_num
CN = args.core_num
IN = args.iter_num

# Fetch data from data files
running_time = []
copy_time    = []
compute_time = []
blocked_time = []
idle_time    = []

for i in range(1, WN + 1):
  print 'Processing the data for worker ' + str(i) + '...'

  data_file =  args.dir + '/' + args.input_base_name + '-' + str(i)

  print 'Opening file ' + data_file
  data = open(data_file, 'r')
  
  running_found = False
  copy_found    = False
  compute_found = False
  blocked_found = False
  idle_found    = False

  for num, line in enumerate(data):
      x =  re.findall('.* (\d+\.\d+)', line)
      if len(x) > 0:
          xnum = float(decimal.Decimal(x[0]))
          if "Running" in line:
              running_time.append(xnum)
              running_found = True
          if "Copy" in line:
              copy_time.append(xnum)
              copy_found = True
          if "Compute" in line:
              compute_time.append(xnum)
              compute_found = True
          elif "Blocked" in line:
              blocked_time.append(xnum)
              blocked_found = True
          elif "Idle" in line:
              idle_time.append(xnum)
              idle_found = True

  if (not running_found or \
      not copy_found    or \
      not compute_found or \
      not blocked_found or \
      not idle_found):
    print 'File ' + data_file + ' does not match the format.'
    exit(0)
  data.close()

# Plot the results in stack bar 
N = WN + 1
P = 5
ind = np.arange(N)
# ind.append(WN + 2)
# ind.append(WN + 3)
width = 0.5

cache_time = [0.709, 0.618, 0.634, 0.6785, 0.599, 0.6396, 0.589, 0.633]
cache_time = map(operator.mul, cache_time, [CN * IN] * WN)
physbam_time = map(operator.sub, compute_time, cache_time)

mpi = [1.1, 0.063, 3.0]

Data = []
# Data.append(map(operator.div, running_time, [CN * IN] * WN))
# Data.append(map(operator.div, compute_time, [CN * IN] * WN))
Data.append(map(operator.div, physbam_time, [CN * IN] * WN))
Data.append(map(operator.div, cache_time, [CN * IN] * WN))
Data.append(map(operator.div, copy_time, [CN * IN] * WN))
Data.append(map(operator.div, blocked_time, [CN * IN] * WN))
Data.append(map(operator.div, idle_time, [CN * IN] * WN))

for i in range(0, len(Data)):
  Data[i].append(np.mean(Data[i]))

Legends = []
# Legends.append('Running')
# Legends.append('Compute')
Legends.append('Nimbus::Compute Raw')
Legends.append('Nimbus::Compute Cache')
Legends.append('Nimbus::Copy')
Legends.append('Nimbus::Blocked')
Legends.append('Nimbus::Idle')
Legends.append('MPI::Compute Raw')
Legends.append('MPI::Pack/Unpack')
Legends.append('MPI::Block')

Colors = []

Colors.append('#8dd3c7')
Colors.append('#fc8d62')
Colors.append('#ffffb3')
Colors.append('#bebada')
Colors.append('w')

Hatch = []
Hatch.append('')
Hatch.append('')
Hatch.append('')
Hatch.append('')
Hatch.append('')



bottom = []
bottom.append([0] * N)

Parts = []

for i in range(0, P):
  p = plt.bar(ind, Data[i], width, color=Colors[i], hatch=Hatch[i], bottom=bottom[i])
  for num, rect in zip(Data[i], p):
    plt.text(rect.get_x() + rect.get_width()/2.,
             rect.get_y()+ rect.get_height()/2.,
             '{:.1f}'.format(num),
             ha='center', va='center')

  bottom.append(map(add, bottom[i], Data[i]))
  Parts.append(p[0])


p = plt.bar(N, mpi[0], width, color=Colors[0], bottom=0)
Parts.append(p[0])
plt.text(p[0].get_x() + p[0].get_width()/2.,
         p[0].get_y() + p[0].get_height()/2.,
         '{:.1f}'.format(mpi[0]),
         ha='center', va='center')

p = plt.bar(N, mpi[1], width, color=Colors[2], bottom=mpi[0])
Parts.append(p[0])
plt.text(p[0].get_x() + p[0].get_width()/2.,
         p[0].get_y() + p[0].get_height()/2.,
         '{:.1f}'.format(mpi[1]),
         ha='center', va='center')

p = plt.bar(N, mpi[2], width, color=Colors[3], bottom=mpi[0] + mpi[1])
Parts.append(p[0])
plt.text(p[0].get_x() + p[0].get_width()/2.,
         p[0].get_y() + p[0].get_height()/2.,
         '{:.1f}'.format(mpi[2]),
         ha='center', va='center')



plt.ylabel('Time (seconds)')

title  = 'PhysBAM Water Simulation Size 256 Cube, 64 uniform partitions, 100 projection iteration\n'
title += '8 c3.2xlarge EC2 workers each with 8 threads, c3.4xlarge controller with 8 assigning threads\n'
# title += 'job done optimization and disabled Nagle\'s algorithm'
# title += 'worker templates activated'
# title += 'all non-sterile jobs templatized'
# title += 'all non-sterile jobs templatized with complex jobs and memoization'

plt.title(title)
xticks = []
for i in range(1, WN + 1):
  xticks.append('W' + str(i))
xticks.append("Nimbus")
xticks.append("MPI")
plt.xticks(np.arange(N + 1) + width/2., xticks )
# plt.yticks(np.arange(0, total[0] , 6))
# plt.legend( (S[0][0], S[1][0], S[2][0], S[3][0], S[4][0]), ('Translator Compute', 'Compute', 'Translator Copy', 'Copy', 'Idle') )


plt.legend(reversed(Parts), reversed(Legends))

plt.ylim(0, 35)

plt.savefig("test.png")
plt.show()


