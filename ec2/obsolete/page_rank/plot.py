#!/usr/bin/python
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.font_manager as fm

# Configure matplot to understand tex grammer.
# rc('text', usetex=True)
rc('font', size=18)

N = 2
P = 4
ind = np.arange(N) * 1.2
width = 0.2
sep   = 0.22

Legends = [
'Spark',
'Nimbus /wo Templates',
'Nimbus /w CT',
'Nimbus /w CT + WT'
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
                 rect.get_y() + rect.get_height() + 1.5,
                 '{:.2f}'.format(num),
                 ha='center', va='center',
                 fontsize='12')
    return p


spark       = [ 50.0, 50.0 ]
template    = [ 6.07, 10.14 ]
template_nb = [ 5.06, 13.33 ]
template_nc = [ 7.06, 19.15 ]

Parts = []

p = plot_bar(spark, ind - 1.5*sep, Colors[0], True)
Parts.append(p[0])

p = plot_bar(template_nc, ind - 0.5*sep, Colors[1], True)
Parts.append(p[0])

p = plot_bar(template_nb, ind + 0.5*sep, Colors[2], True)
Parts.append(p[0])

p = plot_bar(template, ind + 1.5*sep, Colors[3], True)
Parts.append(p[0])

ticks=['400', '800']
plt.xticks(ind+width/2., ticks)
plt.xlabel('Total partitions')
plt.ylabel('Iteration length (seconds)')
axes = plt.gca()
axes.set_ylim([0, 75])

legend_font = fm.FontProperties(size='14')
plt.legend(Parts, Legends,
           ncol=1, loc=2,
           frameon=False, prop=legend_font)

# title  = 'Pagerank for Wikipedia dump'
# plt.title(title, fontsize='15')

# plt.show()
plt.savefig('pr_strong_scale.pdf')
