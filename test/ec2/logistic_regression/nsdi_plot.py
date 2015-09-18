#!/usr/bin/python
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import rc

# Configure matplot to understand tex grammer.
# rc('text', usetex=True)
rc('font', size=24)

N = 4
P = 4
ind = np.arange(N) * 1.2
width = 0.2
sep   = 0.22

Legends = [
'Spark',
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
                 rect.get_y()+ rect.get_height() + 0.1,
                 '{:.2f}'.format(num),
                 ha='center', va='center',
                 fontsize='xx-small')
    return p


spark       = [ 1.90, 2.00, 3.00,  4.00]
no_template = [ 0.27, 0.56, 0.81,  1.19]
c_template  = [ 0.21, 0.38, 0.66,  0.82]
cw_template = [ 0.06, 0.14, 0.20,  0.27]

Parts = []

p = plot_bar(spark, ind - 2*sep, Colors[0], True)
Parts.append(p[0])

p = plot_bar(no_template, ind - sep, Colors[1], True)
Parts.append(p[0])

p = plot_bar(c_template, ind, Colors[2], True)
Parts.append(p[0])

p = plot_bar(cw_template, ind + sep, Colors[3], True)
Parts.append(p[0])

ticks=['20 (160)','40 (320)','60 (480)','80 (640)']
plt.xticks(ind+width/2., ticks)
plt.xlabel('Number of workers (#cores)')
plt.ylabel('Iteration length (seconds)')
plt.ylim([0, 5])
# plt.set_yscale('log')

plt.legend(Parts, Legends,
           ncol=1, loc=1, mode='expand',
           fontsize='small', frameon=False)

title  = 'Logistic Regression, c3.2xlarge worker, c3.4xlarge controller '
title += '5 million 10 dimensional samples per worker'
# plt.title(title, fontsize='small')

plt.show()
# plt.savefig('../figs/weak_scale.pdf')
