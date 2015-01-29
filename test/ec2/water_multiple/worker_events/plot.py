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
blocked_time = []
idle_time = []

for i in range(1, WN + 1):
  print 'Processing the data for worker ' + str(i) + '...'

  data_file =  args.dir + '/' + args.input_base_name + '-' + str(i)

  print 'Opening file ' + data_file
  data = open(data_file, 'r')
  
  running_found = False
  blocked_found = False
  idle_found    = False

  for num, line in enumerate(data):
      x =  re.findall('.* (\d+\.\d+)', line)
      if len(x) > 0:
          xnum = float(decimal.Decimal(x[0]))
          if "Running" in line:
              running_time.append(xnum)
              running_found = True
          elif "Blocked" in line:
              blocked_time.append(xnum)
              blocked_found = True
          elif "Idle" in line:
              idle_time.append(xnum)
              idle_found = True

  if (not running_found or \
      not blocked_found or \
      not idle_found):
    print 'File ' + data_file + ' does not match the format.'
    exit(0)
  data.close()

# Plot the results in stack bar 
N = WN + 2
P = 3
ind = np.arange(N)
# ind.append(WN + 2)
# ind.append(WN + 3)
width = 0.5

Data = []
Data.append(map(operator.div, running_time, [CN * IN] * WN))
Data.append(map(operator.div, blocked_time, [CN * IN] * WN))
Data.append(map(operator.div, idle_time, [CN * IN] * WN))

for i in range(0, len(Data)):
  Data[i].append(np.mean(Data[i]))

Data[0].append(1.1)
Data[1].append(3.0)
Data[2].append(0.0)


Legends = []
Legends.append('Running')
Legends.append('Blocked')
Legends.append('Idle')

Colors = []

Colors.append('#8dd3c7')
# Colors.append('#ffffb3')
Colors.append('#bebada')
Colors.append('w')

Hatch = []
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

plt.ylabel('Time (seconds)')

title  = 'PhysBAM Water Simulation Size 256 Cube, 64 uniform partitions, 100 projection iteration\n'
title += '8 c3.2xlarge EC2 workers each with 8 threads, c3.4xlarge controller with 8 assigning threads\n'
# title += 'job done optimization and disabled Nagle\'s algorithm'
# title += 'worker templates activated'
# title += 'all non-sterile jobs templatized'
title += 'all non-sterile jobs templatized with complex jobs and memoization'

plt.title(title)
xticks = []
for i in range(1, WN + 1):
  xticks.append('W' + str(i))
xticks.append("mean")
xticks.append("PhysBAM")
plt.xticks(ind+width/2., xticks )
# plt.yticks(np.arange(0, total[0] , 6))
# plt.legend( (S[0][0], S[1][0], S[2][0], S[3][0], S[4][0]), ('Translator Compute', 'Compute', 'Translator Copy', 'Copy', 'Idle') )


plt.legend(reversed(Parts), reversed(Legends))

plt.ylim(0, 35)

plt.savefig("test.png")
plt.show()


