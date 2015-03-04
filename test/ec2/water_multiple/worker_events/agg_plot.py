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

# 256 old
# running_time = [5.7,  5.6,  4.2,  4.0 ]
# blocked_time = [2.9,  2.8,  3.3,  2.9 ]
# idle_time    = [21.9, 17.1, 12.8, 15.6]


# 256
running_time = [5.7,  5.6,   5.1,  5.1, 4.9, 3.6]
blocked_time = [2.9,  2.8,   5.7, 10.6, 4.5, 3.5]
idle_time    = [21.9, 17.1, 11.7,  4.4, 5.5, 7.7]

# 512
running_time = [19.3, 14.4]
blocked_time = [16.8, 13.7]
idle_time    = [4.5 , 5.2]

# 256 old 
# physbam = [1.1, 3.0]

# 256
physbam = [1.7, 2.4]

# 512
physbam = [12, 12.4]

# Plot the results in stack bar 
N = len(running_time)
P = 3
ind = np.arange(N)
width = 0.5

Data = []
Data.append(running_time)
Data.append(blocked_time)
Data.append(idle_time)


Legends = []
Legends.append('Nimbus::Running')
Legends.append('Nimbus::Blocked')
Legends.append('Nimbus::Idle')
Legends.append('MPI::Compute')
Legends.append('MPI::Blocked')

Colors = []

Colors.append('#8dd3c7')
# Colors.append('#ffffb3')
Colors.append('#bebada')
Colors.append('w')


bottom = []
bottom.append([0] * N)

Parts = []

for i in range(0, P):
  p = plt.bar(ind, Data[i], width, color=Colors[i], hatch='', bottom=bottom[i])
  for num, rect in zip(Data[i], p):
    plt.text(rect.get_x() + rect.get_width()/2.,
             rect.get_y() + rect.get_height()/2.,
             '{:.1f}'.format(num),
             ha='center', va='center')

  bottom.append(map(add, bottom[i], Data[i]))
  if (i == P - 1):
    for num, rect in zip(bottom[i+1], p):
      plt.text(rect.get_x() + rect.get_width()/2.,
               rect.get_y() + rect.get_height() + .5,
               '{:.1f}'.format(num),
               ha='center', va='center')

  Parts.append(p[0])


p = plt.bar(N, physbam[0], width, color='#ffffb3', bottom=0)
Parts.append(p[0])
plt.text(p[0].get_x() + p[0].get_width()/2.,
         p[0].get_y() + p[0].get_height()/2.,
         '{:.1f}'.format(physbam[0]),
         ha='center', va='center')

p = plt.bar(N, physbam[1], width, color='#fc8d62', bottom=physbam[0])
Parts.append(p[0])
plt.text(p[0].get_x() + p[0].get_width()/2.,
         p[0].get_y() + p[0].get_height()/2.,
         '{:.1f}'.format(physbam[1]),
         ha='center', va='center')

plt.text(p[0].get_x() + p[0].get_width()/2.,
         p[0].get_y() + p[0].get_height() + .5,
         '{:.1f}'.format(physbam[0] + physbam[1]),
         ha='center', va='center')

plt.ylabel('Average Iteration Time (seconds)')
# plt.xlabel('Different Simulation Setups')

# title  = 'PhysBAM Water Simulation Size 256 Cube, 64 uniform partitions, 100 projection iteration\n'
title  = 'PhysBAM Water Simulation Size 512 Cube, 64 uniform partitions, 100 projection iteration\n'
title += '8 c3.2xlarge EC2 workers each with 8 threads, c3.4xlarge controller with 8 assigning threads\n'

plt.title(title)

# 256
xticks = ['No Template', 'Worker Template \n (Only Projection)', 'Complete Worker Template \n + Complex Jobs', 'Binding Template', 'No Explicit \n Inter-Worker Before Set', 'Removed Worker Logging', 'PhysBAM']

# 512
xticks = ['Temlates and Memoization', 'Removed Worker Logging', 'PhysBAM']

plt.xticks(np.arange(N + 1) + width/2., xticks )
# plt.yticks(np.arange(0, total[0] , 6))
# plt.legend( (S[0][0], S[1][0], S[2][0], S[3][0], S[4][0]), ('Translator Compute', 'Compute', 'Translator Copy', 'Copy', 'Idle') )


plt.legend(reversed(Parts), reversed(Legends))

# 256
plt.ylim(0, 35)

# 512
plt.ylim(0, 45)

plt.savefig("test.png")
plt.show()


