#!/usr/bin/python
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import rc

# Configure matplot to understand tex grammer.
# rc('text', usetex=True)
rc('font', size=24)

N = 5
P = 4
ind = np.arange(N) * 1.2
width = 0.2
sep   = 0.22

Legends = [
'Spark',
'Nimbus /wo Template',
'Nimbus /w Template (Including Learning)',
'Nimbus /w Template (Excluding Learning)'
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
                 rect.get_y()+ rect.get_height() + 0.3,
                 '{:.2f}'.format(num),
                 ha='center', va='center',
                 fontsize='xx-small')
    return p


spark       = [ 7.30, 7.50, 7.40, 10.50, 22.20]
no_template = [ 0.20, 0.65, 1.41,  2.92,  7.76]
template_il = [ 0.14, 0.24, 0.45,  0.93,  2.34]
template_xl = [ 0.14, 0.15, 0.24,  0.45,  1.11]

Parts = []

p = plot_bar(spark, ind - 2*sep, Colors[0], True)
Parts.append(p[0])

p = plot_bar(no_template, ind - sep, Colors[1], True)
Parts.append(p[0])

p = plot_bar(template_il, ind, Colors[2], True)
Parts.append(p[0])

p = plot_bar(template_xl, ind + sep, Colors[3], True)
Parts.append(p[0])

ticks=['4','20','40','80','200']
plt.xticks(ind+width/2., ticks)
plt.xlabel('Partition number per core (200 cores)')
plt.ylabel('Iteration length (seconds)')
# plt.set_yscale('log')

plt.legend(Parts, Legends,
           ncol=1, loc=1, mode='expand',
           fontsize='small', frameon=False)

title  = 'Logistic Regression, 25 workers, 1 controller, m3.2xlarge instances '
title += '400 million 10 dimensional samples'
plt.title(title, fontsize='small')

plt.show()
# plt.savefig('../figs/weak_scale.pdf')
