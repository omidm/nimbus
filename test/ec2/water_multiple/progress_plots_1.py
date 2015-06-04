#!/usr/bin/python
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import rc

# Configure matplot to understand tex grammer.
# rc('text', usetex=True)
rc('font', size=24)

compute_time = [10.91,  11.0, 11.07, 11.05, 11.10]
cache_time   = [ 2.17,  2.23,  2.23,  2.19,  2.24]
block_time   = [13.28, 13.27, 13.10, 12.66, 12.47]
idle_time    = [ 4.33,  4.03,  3.84,  3.50,  3.49]


# smoothe compute and cache time among all of them.
compute_time = [np.mean(compute_time)] * len(compute_time)
cahe_time = [np.mean(cache_time)] * len(cache_time)


Data = []
Data.append(compute_time)
Data.append(cache_time)
Data.append(block_time)
Data.append(idle_time)

RefData = [[12.31], [.25], [13.00]]


N = len(Data[0])
ind = np.arange(N)
width = 0.4




Legends = [
'PhysBAM:Compute',
'PhysBAM:Serialize/Deserialize',
'PhysBAM:Blocked',
'Nimbus:Compute',
'Nimbus:Copy/Cache',
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



ytop = 4
labels = ['Cache+\nTemplate',
          'Exclude\nInit Idle',
          'No Job Done\nFlooding',
          'Priority Queue\nFor Copy Jobs',
          'Multi-Threaded\nAfter Set Clean Up',
          'VI', 'VII', 'VIII']
ticks= labels[0 : N] + ['PhysBAM']
plt.xticks(np.append(ind, [N]) +width/2., ticks)
plt.yticks(np.linspace(0,ytop*10,ytop+1))
title  = 'PhysBAM Water Simulation Size 512 Cube, 64 uniform partitions, 100 projection iteration\n'
title += '8 c3.2xlarge EC2 workers each with 8 threads, c3.4xlarge controller with 8 assigning threads\n'
plt.title(title)
plt.ylabel('Time(s)')

plt.xlim(-0.3,N+.7)
plt.ylim(0, ytop*10)

plt.grid(axis='y')

plt.legend(reversed(Parts), reversed(Legends),
           ncol=2, loc='upper right',
           fontsize='medium', frameon=False)


plt.show()
