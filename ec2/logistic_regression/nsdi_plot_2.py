#!/usr/bin/python
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import rc

# Configure matplot to understand tex grammer.
# rc('text', usetex=True)
rc('font', size=24)

N = 3
P = 4
ind = np.arange(N) * 1.2
width = 0.25
sep   = 0.27

Legends = [
'Spark',
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
                 rect.get_y()+ rect.get_height() + 0.1,
                 '{:.2f}'.format(num),
                 ha='center', va='center',
                 fontsize='small')
    return p


spark       = [ 1.54, 1.56, 1.64]
no_template = [ 0.24, 0.54, 1.12]
c_template  = [ 0.19, 0.40, 0.84]
cw_template = [ 0.06, 0.12, 0.22]

Parts = []

p = plot_bar(spark, ind - 2*sep, RColors[0], True)
Parts.append(p[0])

p = plot_bar(no_template, ind - sep, Colors[1], True)
Parts.append(p[0])

p = plot_bar(c_template, ind, Colors[2], True)
Parts.append(p[0])

p = plot_bar(cw_template, ind + sep, Colors[3], True)
Parts.append(p[0])

ticks=['(25, 2000, 100)','(25, 4000, 100)', '(25, 8000, 100)']
plt.xticks(ind, ticks, fontsize='medium')
plt.xlabel('(#workers, #partitions, #samples in millions)', fontsize='large')
plt.ylabel('Iteration length (seconds)', fontsize='large')
plt.ylim([0, 3])
plt.yticks(np.linspace(0,3,4))
plt.xlim([-1 + width, N + width])
# plt.set_yscale('log')

plt.legend(Parts, Legends,
           ncol=1, loc=1, mode='expand',
           fontsize='large', frameon=False)

title  = 'Logistic Regression, c3.2xlarge worker, c3.4xlarge controller '
title += '5 million 10 dimensional samples per worker'
# plt.title(title, fontsize='small')

plt.show()
# plt.savefig('../figs/weak_scale.pdf')
