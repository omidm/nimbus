#!/usr/bin/python
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import rc

# Configure matplot to understand tex grammer.
# rc('text', usetex=True)
rc('font', size=24)

run_time   = [13.09, 13.04]
block_time = [24.58, 20.76]
idle_time  = [11.47, 11.62]


# smoothe run and cache time among all of them.
run_time = [np.mean(run_time)] * len(run_time)


Data = []
Data.append(run_time)
Data.append(block_time)
Data.append(idle_time)

RefData = [[12.16+.26], [19.27]]


N = len(Data[0])
ind = np.arange(N)
width = 0.4




Legends = [
'PhysBAM:Compute',
'PhysBAM:Blocked',
'Nimbus:Compute',
'Nimbus:Blocked',
'Nimbus:Idle',
]

n_colors = 3
color_map = cmx.ScalarMappable(
        colors.Normalize(vmin=0, vmax=n_colors+1), cmap='Greens')
Colors = [color_map.to_rgba(i) for i in range(n_colors, -1, -1)]

color_map = cmx.ScalarMappable(
        colors.Normalize(vmin=0, vmax=n_colors+1), cmap='Reds')
SubColors = [color_map.to_rgba(i) for i in range(n_colors, 0, -1)]

def plot_bar(bar_data, ind, bottom, color, hatch, add_sum):
    # add_sum = False
    p = plt.bar(ind, bar_data, width,
                color=color, hatch=hatch, bottom=bottom)
    i = 0
    for num, rect in zip(bar_data, p):
      plt.text(rect.get_x() + rect.get_width()/2.,
               rect.get_y()+ rect.get_height()/2.,
               '{:.1f}'.format(num),
               ha='center', va='center',
               fontsize='xx-small')
      if add_sum:
          plt.text(rect.get_x() + rect.get_width()/2,
                   rect.get_y()+ rect.get_height() + 0.5,
                   '{:.1f}'.format(num+bottom[i]),
                   ha='center', va='center',
                   fontsize='xx-small')
          i += 1
    for i in range(0, len(bottom)):
        bottom[i] += bar_data[i]
    return p



Parts = []

bottom = [0] * 1
for i in range(0, len(RefData)):
    p = plot_bar(RefData[i], [N], bottom, SubColors[i], '', i==len(RefData)-1)
    Parts.append(p[0])

bottom = [0] * N
for i in range(0, len(Data)):
    p = plot_bar(Data[i], ind, bottom, Colors[i], '', i==len(Data)-1)
    Parts.append(p[0])



ytop = 6
labels = ['After All Previous\nBlock Time Optimizations',
          'TCP Buffer Tuning+\nIncreased Receive Event Quota',
          'VII', 'VIII']
ticks= labels[0 : N] + ['PhysBAM\nExcluded Init/Write']
plt.xticks(np.append(ind, [N]) +width/2., ticks)
plt.yticks(np.linspace(0,ytop*10,ytop+1))
title  = 'PhysBAM Water Simulation Size 1024 Cube, 512 uniform partitions, 100 projection iteration\n'
title += '64 c3.2xlarge EC2 workers each with 8 threads, c3.4xlarge controller with 8 assigning threads\n'
plt.title(title, fontsize="medium")
plt.ylabel('Time(s)')

plt.xlim(-0.3,N+.7)
plt.ylim(0, ytop*10)

plt.grid(axis='y')

plt.legend(reversed(Parts), reversed(Legends),
           ncol=2, loc='upper right',
           fontsize='medium', frameon=False)


plt.show()
