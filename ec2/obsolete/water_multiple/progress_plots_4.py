#!/usr/bin/python
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import rc

# Configure matplot to understand tex grammer.
# rc('text', usetex=True)
rc('font', size=24)

run_time = [13.14, 13.29, 13.18]

run_time = [np.mean(run_time)] * len(run_time)

Data = [
# run time
[
run_time,              # avg(run)
],

# synch time
[
[19.89, 22.42, 21.34], # avg(block)
[ 3.30,  2.11,  2.12], # avg(idle) - spawner_worker(idle)
],

# idle time
[
[ 4.65,  3.11,  2.91], # controller(AC+DMQ+TI)
[ 1.24,  0.85,  0.83], # parent_exec
[ 0.99,  0.52,  0.48]  # spawner_worker(idle) - controller(AC+DMQ+TI) - parent_exec
]

]


RefData = [[12.16+.26], [19.27]]



N = len(Data[0][0])
ind = np.arange(N)
width = 0.4




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





def plot_bar(bar_data, ind, bottom, color, hatch, add_sum):
    # add_sum = False
    p = plt.bar(ind, bar_data, width,
                color=color, hatch=hatch, bottom=bottom)

    i = 0
    for num, rect in zip(bar_data, p):
      if (num > 0.5):
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
    for j in range(0, len(Data[i])):
      p = plot_bar(Data[i][j], ind, bottom, Colors[i], Hatch[j], (j==len(Data[i])-1) and (i==len(Data)-1))
      Parts.append(p[0])









ytop = 6
labels = ['Base\nNo Batching',
          'Projection Batching\nSize 10',
          'Remove Pjoection\nBottleneck Job',
          'VIII']
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
           ncol=3, loc='upper right',
           fontsize='medium', frameon=False)


plt.show()











# ytop = 5
# ticks=['512 (64)','1024 (512)','worker 3','worker 4','worker 5', 'worker 6', 'worker 7', 'worker 8', 'Average']
# plt.xticks(ind+width/2., ticks)
# plt.yticks(np.linspace(0,ytop*10,ytop+1))
# plt.xlabel('Simulation Size (partition #)')
# plt.ylabel('Time(s)')
# 
# margin = 0.2
# plt.xlim(-0.2-margin,N-0.5+.6+margin)
# # plt.ylim(0,35)
# plt.legend(reversed(Parts), reversed(Legends),
#            ncol=3, loc=3, bbox_to_anchor=(0.,1.02,1.,.102),
#            borderaxespad=0., mode='expand',
#            fontsize='small', frameon=False)
# 
# plt.gca().spines['top'].set_visible(False)
# plt.gca().spines['right'].set_visible(False)
# plt.grid(axis='y')
# plt.gca().get_xaxis().set_tick_params(top='off')
# plt.tight_layout(pad=0.1, rect=(0,0,1,0.8))
# 
# plt.show()
# # plt.savefig('../figs/weak_scale.pdf')




