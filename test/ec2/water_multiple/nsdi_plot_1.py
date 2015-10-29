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
'MPI',
'Nimbus /wo Templates',
'Nimbus /w  CT',
'Nimbus /w  CT + WT'
]

n_colors = P
color_map = cmx.ScalarMappable(
        colors.Normalize(vmin=0, vmax=n_colors+1), cmap='Greens')
Colors = [color_map.to_rgba(i) for i in range(n_colors, -1, -1)]

color_map = cmx.ScalarMappable(
        colors.Normalize(vmin=0, vmax=n_colors+1), cmap='Reds')
RColors = [color_map.to_rgba(i) for i in range(n_colors, -1, -1)]


def plot_bar(bar_data, ind, color, add_sum):
    p = plt.bar(ind, bar_data, width, color=color)
    if add_sum:
      for num, rect in zip(bar_data, p):
        plt.text(rect.get_x() + rect.get_width()/2,
                 rect.get_y()+ rect.get_height() + 3,
                 '{:.1f}'.format(num),
                 ha='center', va='center',
                 fontsize='small')
    return p


physbam     = [24.70,  31.70]
no_template = [35.21, 196.78]
c_template  = [27.75, 173.28]
cw_template = [25.32,  36.50]
# cw_template = [25.32,  37.90] no batching for 1024

Parts = []

p = plot_bar(physbam, ind - 2*sep, RColors[0], True)
Parts.append(p[0])

p = plot_bar(no_template, ind - sep, Colors[1], True)
Parts.append(p[0])

p = plot_bar(c_template, ind, Colors[2], True)
Parts.append(p[0])

p = plot_bar(cw_template, ind + sep, Colors[3], True)
Parts.append(p[0])

ticks=['(8, 64, 512$^3$)','(64, 512, 1024$^3$)']
plt.xticks(ind, ticks, fontsize='medium')
plt.xlabel('(#workers, #partitions, #cells)', fontsize='large')
plt.ylabel('Iteration length (seconds)', fontsize='large')
plt.ylim([0, 210])
# plt.set_yscale('log')

plt.legend(Parts, Legends,
           ncol=1, loc=1, mode='expand',
           fontsize='large', frameon=False)

title  = 'Water Simulation, c3.2xlarge worker, c3.4xlarge controller '
# plt.title(title, fontsize='small')

plt.show()
# plt.savefig('../figs/weak_scale.pdf')
