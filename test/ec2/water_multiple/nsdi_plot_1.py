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
P = 4
ind = np.arange(N) * 1.2
width = 0.2
sep   = 0.22

Legends = [
'PhysBAM',
'Nimbus /wo Templates',
'Nimbus /w  Controller Template',
'Nimbus /w  Controller + Worker Template'
]

n_colors = P
color_map = cmx.ScalarMappable(
        colors.Normalize(vmin=0, vmax=n_colors+1), cmap='Greens')
Colors = [color_map.to_rgba(i) for i in range(n_colors, -1, -1)]

def plot_bar(bar_data, ind, color, add_sum):
    p = plt.bar(ind, bar_data, width, color=color)
    if add_sum:
      for num, rect in zip(bar_data, p):
        plt.text(rect.get_x() + rect.get_width()/2,
                 rect.get_y()+ rect.get_height() + 3,
                 '{:.1f}'.format(num),
                 ha='center', va='center',
                 fontsize='xx-small')
    return p


physbam     = [24.70,  31.70]
no_template = [35.21, 196.78]
c_template  = [27.75, 173.28]
cw_template = [25.32,  37.90]

Parts = []

p = plot_bar(physbam, ind - 2*sep, Colors[0], True)
Parts.append(p[0])

p = plot_bar(no_template, ind - sep, Colors[1], True)
Parts.append(p[0])

p = plot_bar(c_template, ind, Colors[2], True)
Parts.append(p[0])

p = plot_bar(cw_template, ind + sep, Colors[3], True)
Parts.append(p[0])

ticks=['8 (64)','64 (512)']
plt.xticks(ind+width/2., ticks)
plt.xlabel('Number of workers (#cores)')
plt.ylabel('Iteration length (seconds)')
plt.ylim([0, 210])
# plt.set_yscale('log')

plt.legend(Parts, Legends,
           ncol=1, loc=1, mode='expand',
           fontsize='small', frameon=False)

title  = 'Water Simulation, c3.2xlarge worker, c3.4xlarge controller '
# plt.title(title, fontsize='small')

plt.show()
# plt.savefig('../figs/weak_scale.pdf')
