import decimal
import argparse
import os
import os.path
import re

import numpy as np
import matplotlib.pyplot as plt
from operator import add
import operator

def parse_line(line):
    items = line.split()
    return items[3], float(items[5])

stats = list()
for rank in range(1, 9):
    f = open('{}_time_per_thread.txt'.format(rank), 'r')
    stat = dict()
    for line in f:
        key, time = parse_line(line)
        if key in stat:
            stat[key] += time / 8 / 55
        else:
            stat[key] = time / 8 / 55
    #print rank, stat
    #print rank, stat['kExecuteCopyJob'] + stat['kExecuteComputationJob']
    stats.append(stat)

print eval(str(stats))
# Worker number.
WN = 8
# Number of bars.
N = WN + 1
ind = np.arange(N)
width = 0.3

def StatsToData(stats):
    Data = []
    Data.append(
            [item['kExecuteComputationJob'] - item['kAssemblingCache']
             for item in stats])
    Data.append(
            [item['kAssemblingCache']
             for item in stats])
    Data.append(
            [item['kExecuteCopyJob']
             for item in stats])
    Data.append(
            [item['kTotal']
             - item['kExecuteCopyJob'] - item['kExecuteComputationJob']
             for item in stats])
    # Append mean
    for i in range(0, len(Data)):
      Data[i].append(np.mean(Data[i]))
    return Data

Data = StatsToData(stats)
RefData = list()
for item in Data:
    RefData.append(item[:])

#RefData[0] = [1.6,  .4, 1.4, 1.6, 0.5, 1.5,  .6,  .6, 1.0, 1.1]
#RefData[1] = [ .7,  .6,  .6,  .7,  .6,  .6,  .6,  .6,  .6,  .1]
#RefData[2] = [3.1, 1.9, 2.8, 2.8, 1.9, 2.5, 1.7, 1.7, 2.3,   0]
#x = [3.7, 2.9, 2.8, 3.0, 3.0, 3.7, 2.9, 3.2, 3.2, 3.0]
#y = [13.9,17.2,15.2,14.9,16.8,14.6,17.2,16.8,15.8,0]
#RefData[3] = [sum(i) for i in zip(x,y)]

# comp, pack+unpack, block
PhysbamData = [[1.65631185292, 0.0307991847995+0.0311240407941, 2.27145996632],
               [1.6786182701,  0.0330878960768+0.0320284208514, 2.45848018238]]

Legends = []
Legends.append('Nimbus::Compute Raw')
Legends.append('Nimbus::Compute Cache')
Legends.append('Nimbus::Copy')
#Legends.append('Nimbus::Blocked')
Legends.append('Nimbus::Idle+Blocked')
Legends.append('MPI::Compute Raw')
Legends.append('MPI::Pack/Unpack')
Legends.append('MPI::Block')

Colors = []
Colors.append('#8dd3c7')
Colors.append('#fdc086')
Colors.append('#ffffb3')
#Colors.append('#bebada')
Colors.append('w')

Hatch = []
Hatch.append('')
Hatch.append('')
Hatch.append('')
Hatch.append('')
Hatch.append('')

def plot_bar(bar_data, ind, bottom, color, hatch, left):
    if left:
        p = plt.bar(ind - width/2, bar_data, width,
                    color=color, hatch=hatch, bottom=bottom)
    else:
        p = plt.bar(ind + width/2, bar_data, width,
                    color=color, hatch=hatch, bottom=bottom)
    for num, rect in zip(bar_data, p):
      plt.text(rect.get_x() + rect.get_width()/2.,
               rect.get_y()+ rect.get_height()/2.,
               '{:.1f}'.format(num),
               ha='center', va='center')
    for i in range(0, len(bottom)):
        bottom[i] += bar_data[i]
    return p



bottom = [0] * N
for i in range(0, len(RefData)):
    plot_bar(RefData[i], ind, bottom, Colors[i], Hatch[i], False)

bottom = [0]
NewColors = [Colors[0], Colors[1], Colors[3]]
for i in range(0, len(PhysbamData[0])):
    plot_bar([PhysbamData[0][i]],
             np.array([WN+1]),
             bottom, NewColors[i], Hatch[i], False)

Parts = []
bottom = [0] * N
for i in range(0, len(Data)):
    p = plot_bar(Data[i], ind, bottom, Colors[i], Hatch[i], True)
    Parts.append(p[0])

bottom = [0]
NewColors = [Colors[0], Colors[1], Colors[3]]
for i in range(0, len(PhysbamData[0])):
    p = plot_bar([PhysbamData[1][i]],
                 np.array([WN+1]),
                 bottom, NewColors[i], Hatch[i], True)
    Parts.append(p[0])

title  = 'PhysBAM Water Simulation Size 256 Cube, 64 uniform partitions, 100 projection iteration\n'
title += '8 c3.2xlarge EC2 workers each with 8 threads, c3.4xlarge controller with 8 assigning threads\n'
title += 'worker templates disabled'
plt.title(title)

xticks = []
for i in range(1, WN + 1):
  xticks.append('W' + str(i) + ' (old/new)')
xticks.append("mean (old/new)")
xticks.append("PhysBAM\n(worst in 6 samples/average)")
ind = np.arange(N+1)
plt.xticks(ind+width/2., xticks )
plt.ylabel('Time (seconds)')
plt.xlim(-0.5, N+1)
plt.ylim(0, 40)

plt.legend(reversed(Parts), reversed(Legends))

plt.show()
