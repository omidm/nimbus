#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def plot_nimbus(data, name, bottom, color):
    rects = ax.bar(ind, data, width, bottom=bottom, color=color, label=name)
    for i in range(0, len(bottom)):
        bottom[i] = bottom[i] + data[i]
    for num, rect in zip(data, rects):
        ax.text(rect.get_x()+rect.get_width()/2.,
                rect.get_y()+ rect.get_height()/2.,
                '{:.1f}'.format(num), ha='center', va='center')

def plot_physbam(data, name, bottom, color):
    rects = ax.bar(ind + width, data, width, bottom=bottom, color=color, label=name)
    for i in range(0, len(bottom)):
        bottom[i] = bottom[i] + data[i]
    for num, rect in zip(data, rects):
        ax.text(rect.get_x()+rect.get_width()/2.,
                rect.get_y()+ rect.get_height()/2.,
                '{:.1f}'.format(num), ha='center', va='center')

N = 6
f = open('summary.txt')
lines = f.readlines()
# comp_time, non_comp_time, ready_time, blocked_time, idle_time.
nimbus_data = []
for rank in range(0, N*2, 2):
    nimbus_data.append(eval(lines[rank]))
f.close()
f = open('../../physbam/temp_analyze/summary.txt', 'r')
lines = f.readlines()
physbam_data = []
for rank in range(0, N):
    physbam_data.append([float(lines[rank].split()[0]),
                         float(lines[rank].split()[1])])
print nimbus_data, physbam_data

ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars

fig, ax = plt.subplots()
bottom = [0] * N

plot_nimbus([item[0] + item[2] for item in nimbus_data], 'nimbus:compute_time',
            bottom, 'g')
plot_nimbus([item[3] for item in nimbus_data], 'nimbus:block_time', bottom, 'r')
plot_nimbus([item[1] for item in nimbus_data], 'nimbus:copy_time', bottom, 'b')
plot_nimbus([item[4] for item in nimbus_data], 'nimbus:idle_time', bottom, 'w')

bottom = [0] * N
plot_physbam([item[0] for item in physbam_data], 'physbam:compute_time', bottom,
             'g')
plot_physbam([item[1] for item in physbam_data], 'physbam:block_time', bottom,
             'r')

ax.set_xlabel('Scale')
ax.set_ylabel('Time/Iteration(Second)')
ax.set_title('Settings: 8 workers, c3.2xlarge machine, 8 threads each worker.\n'\
             '3 frames, and projeciton iteration threadhold is 100.')
ax.set_xticks(ind+width)
scale = [256, 344 ,400, 440, 480, 512]
ax.set_xticklabels([str(i) for i in scale])

ax.legend()

plt.show()
