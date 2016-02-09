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
    type=int,
    help="number of workers to process")

args = parser.parse_args()

N = args.worker_num
P = 8

# Fetch data from data files
compute_non_sterile_bare = []
compute_non_sterile_trans = []
compute_sterile_bare = []
compute_sterile_trans = []
copy_bare = []
copy_trans = []
projection_bare = []
projection_trans = []
idle = []
total = []

for i in range(1, N + 1):
  print 'Processing the data for worker ' + str(i) + '...'

  data_file =  args.dir + '/' + args.input_base_name + '-' + str(i)

  print 'Opening file ' + data_file
  data = open(data_file, 'r')
  
  cnb_found = False
  cnt_found = False
  csb_found = False
  cst_found = False
  cob_found = False
  cot_found = False
  prb_found = False
  prt_found = False
  i_found   = False
  t_found   = False

  for num, line in enumerate(data):
      x =  re.findall('(\d+\.\d+$|\d+e-\d+$|\d+$)', line)
      if len(x) > 0:
          xnum = decimal.Decimal(x[0])
          if "Application" in line:
              total.append(xnum)
              t_found = True
          if "Decomposed Compute Non-Sterile Bare" in line:
              compute_non_sterile_bare.append(xnum)
              cnb_found = True
          if "Decomposed Compute Non-Sterile Translator" in line:
              compute_non_sterile_trans.append(xnum)
              cnt_found = True
          if "Decomposed Compute Sterile Bare" in line:
              compute_sterile_bare.append(xnum)
              csb_found = True
          if "Decomposed Compute Sterile Translator" in line:
              compute_sterile_trans.append(xnum)
              cst_found = True
          if "Decomposed Copy Bare" in line:
              copy_bare.append(xnum)
              cob_found = True
          if "Decomposed Copy Translator" in line:
              copy_trans.append(xnum)
              cot_found = True
          if "Decomposed Projection Bare" in line:
              projection_bare.append(xnum)
              prb_found = True
          if "Decomposed Projection Translator" in line:
              projection_trans.append(xnum)
              prt_found = True
          if "Decomposed Idle" in line:
              idle.append(xnum)
              i_found = True

  if (not cnb_found or \
      not cnt_found or \
      not csb_found or \
      not cst_found or \
      not cob_found or \
      not cot_found or \
      not prb_found or \
      not prt_found or \
      not i_found or \
      not t_found):
    print 'File ' + data_file + ' does not match the format.'
    exit(0)
  data.close()

# Plot the results in stack bar 
ind = np.arange(N)
width = 0.35

Data = []
Data.append(map(operator.sub, compute_sterile_bare, projection_bare))
Data.append(projection_bare)
Data.append(compute_sterile_trans)
Data.append(copy_bare)
Data.append(copy_trans)
Data.append(compute_non_sterile_bare)
Data.append(compute_non_sterile_trans)
Data.append(idle)


print compute_sterile_bare

Legends = []
Legends.append('Compute Sterile')
Legends.append('Compute Sterile Projection')
Legends.append('Compute Sterile Translator')
Legends.append('Copy')
Legends.append('Copy Translator')
Legends.append('Compute Non-Sterile')
Legends.append('Compute Non-Sterile Translator')
Legends.append('Idle')

Colors = []

Colors.append('#8dd3c7')
Colors.append('#8dd3c7')
Colors.append('#8dd3c7')
Colors.append('#ffffb3')
Colors.append('#ffffb3')
Colors.append('#bebada')
Colors.append('#bebada')
Colors.append('w')

# Colors.append('#4040ff')
# Colors.append('#4040ff')
# Colors.append('#00bfff')
# Colors.append('#00bfff')
# Colors.append('#afeeee')
# Colors.append('#afeeee')
# Colors.append('w')


Hatch = []
Hatch.append('')
Hatch.append('\\')
Hatch.append('/')
Hatch.append('')
Hatch.append('/')
Hatch.append('')
Hatch.append('/')
Hatch.append('')



bottom = []
bottom.append([0] * N)

Parts = []

for i in range(0, P):
  p = plt.bar(ind, Data[i], width, color=Colors[i], hatch=Hatch[i], bottom=bottom[i])
  print Data[i]
  bottom.append(map(add, bottom[i], Data[i]))
  Parts.append(p[0])

plt.ylabel('Time (seconds)')
plt.title('PhysBAM Water Simulation Size 256 Cube with 8 Nimbus Worker Nodes')
xticks = []
for i in range(1, N + 1):
  xticks.append('W' + str(i))
plt.xticks(ind+width/2., xticks )
# plt.yticks(np.arange(0, total[0] , 6))
# plt.legend( (S[0][0], S[1][0], S[2][0], S[3][0], S[4][0]), ('Translator Compute', 'Compute', 'Translator Copy', 'Copy', 'Idle') )

# annotate the graph with physbam-mpi completion time
physbam = 62 + 118 + 152 + 192 + 229 + 272 + 321 + 402 + 639 + 539 + 598
hline = plt.axhline(physbam, color='#fb8072', ls='--', linewidth=4)
# corner = (ind[N - 2] + ind[N - 3]) / 2.0
# plt.annotate('PhysBAM-MPI\nCompletion Time',
#                      xy=(corner, physbam * 1.6),
#                      xytext=(corner, physbam * 1.6),
#                      color='#fb8072')
Parts.append(hline)
Legends.append('PhysBAM-MPI Completion Time')


plt.legend(reversed(Parts), reversed(Legends))

plt.savefig("test.png")
plt.show()


