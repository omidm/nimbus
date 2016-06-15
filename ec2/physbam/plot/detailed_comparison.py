#!/usr/bin/python
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import rc

# Configure matplot to understand tex grammer.
# rc('text', usetex=True)
rc('font', size=24)

N = 2
ind = np.arange(N - 1)
ind = np.append(ind, [N - 0.5])
width = 0.4



#### Data: no batching for projection (June, 2015)
# Data = [
# # run time
# [
# [12.59               , 13.14             ], # avg(run)
# ],
# 
# # synch time
# [
# [11.69               , 19.89             ], # avg(block)
# [2.56 - 1.22         , 10.18 - 6.88      ], # avg(idle) - min_worker(idle)
# ],
# 
# # idle time
# [
# [0.68               ,  4.65              ], # controller(AC+DMQ+TI)
# [0.14                ,  1.24             ], # parent_exec
# [1.22 - 0.68 - 0.14  , 6.88 - 4.65 - 1.24]  # min_worker(idle) - controller(AC+DMQ+TI) - parent_exec
# ]
# 
# ]


# ### Data: only for 1024, batching projection size 10, and no projection bottleneck (June 2015)
# Data = [
# # run time
# [
# [12.59,13.18], # avg(run)
# ],
# 
# # synch time
# [
# [11.69, 21.34], # avg(block)
# [ 1.34,  2.12], # avg(idle) - min_worker(idle)
# ],
# 
# # idle time
# [
# [ 0.68, 2.91], # controller(AC+DMQ+TI)
# [ 0.14, 0.83], # parent_exec
# [ 0.40, 0.48]  # min_worker(idle) - controller(AC+DMQ+TI) - parent_exec
# ]
# 
# ]


### Data: only for 1024, batching projection size 10, and no projection bottleneck (September 2015)
Data = [
# run time
[
[11.90,12.20], # avg(run)
],

# synch time
[
[10.99, 17.84], # avg(block)
[ 1.10,  1.58], # avg(idle) - min_worker(idle)
],

# idle time
[
[ 0.67, 3.38], # controller(AC+DMQ+TI)
[ 0.14, 0.83], # parent_exec
[ 0.52, 0.71]  # min_worker(idle) - controller(AC+DMQ+TI) - parent_exec
]

]


RefData = [
# run time
[
[11.85, 12.42], # avg(run)
],

# synch time
[
[12.88, 19.27], # avg(block)
]

]



Legends = [
'PhysBAM:Compute',
'PhysBAM:Synch',
'Nimbus:Compute',
'Nimbus:Synch(Blocked)',
'Nimbus:Synch(End of Batch)',
'Nimbus:Idle(Controller Critical)',
'Nimbus:Idle(Parent Execution)',
'Nimbus:Idle(Messaging)',
]

n_colors = 3
color_map = cmx.ScalarMappable(
        colors.Normalize(vmin=0, vmax=n_colors+1), cmap='Greens')
Colors = [color_map.to_rgba(i) for i in range(n_colors, -1, -1)]

color_map = cmx.ScalarMappable(
        colors.Normalize(vmin=0, vmax=n_colors+1), cmap='Reds')
SubColors = [color_map.to_rgba(i) for i in range(n_colors, 0, -1)]

Hatch = ['', '\\\\', '//']





def plot_bar(bar_data, ind, bottom, color, hatch, left, add_sum):
    # add_sum = False
    if left:
        p = plt.bar(ind - width/2, bar_data, width,
                    color=color, hatch=hatch, bottom=bottom)
    else:
        p = plt.bar(ind + width/2, bar_data, width,
                    color=color, hatch=hatch, bottom=bottom)
    i = 0
    for num, rect in zip(bar_data, p):
      if (num > 0.9):
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
bottom = [0] * N
for i in range(0, len(RefData)):
    for j in range(0, len(RefData[i])):
      p = plot_bar(RefData[i][j], ind, bottom, SubColors[i], Hatch[j], False, (j==len(RefData[i])-1) and (i==len(RefData)-1))
      Parts.append(p[0])


bottom = [0] * N
for i in range(0, len(Data)):
    for j in range(0, len(Data[i])):
      p = plot_bar(Data[i][j], ind, bottom, Colors[i], Hatch[j], True, (j==len(Data[i])-1) and (i==len(Data)-1))
      Parts.append(p[0])







ytop = 5
ticks=['512 (64)','1024 (512)']
plt.xticks(ind+width/2., ticks)
plt.yticks(np.linspace(0,ytop*10,ytop+1))
plt.xlabel('Simulation Size (partition #)')
plt.ylabel('Time(s)')

margin = 0.2
plt.xlim(-0.2-margin,N-0.5+.6+margin)
# plt.ylim(0,35)
plt.legend(reversed(Parts), reversed(Legends),
           ncol=3, loc=3, bbox_to_anchor=(0.,1.02,1.,.102),
           borderaxespad=0., mode='expand',
           fontsize='small', frameon=False)

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.grid(axis='y')
plt.gca().get_xaxis().set_tick_params(top='off')
plt.tight_layout(pad=0.1, rect=(0,0,1,0.8))

plt.show()
# plt.savefig('../figs/weak_scale.pdf')




