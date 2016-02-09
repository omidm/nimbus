#!/usr/bin/python
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import rc

# Configure matplot to understand tex grammer.
# rc('text', usetex=True)
rc('font', size=24)

N = 9
ind = np.arange(N - 1)
ind = np.append(ind, [N - 0.5])
width = 0.4

Data = [
[16.16, 15.36, 5.92, 6.22, 16.09, 15.63,  6.03,  5.87, 10.91],
[ 2.77,  2.48, 1.61, 1.60,  2.92,  2.56,  1.81,  1.61,  2.17],
[9.51, 9.64, 17.16, 16.83, 9.02,   9.55, 17.12, 17.39, 13.28],
[2.47, 3.32,  5.94,  5.97, 2.62,  2.88,  5.67,   5.76,  4.33]
# [10.86, 10.51, 10.42, 10.41, 10.20],
# [ 1.72,  1.53, 1.92,   2.18,  2.43],
# [ 5.69,  7.98, 9.75,  11.53, 11.75],
# [ 1.02,  0.99, 1.61,   2.09,  2.59],
]
RefData = [
[17.58+0.28, 18.53+0.32, 6.70+0.17, 7.10+0.19, 17.52+0.32, 17.64+0.33, 6.66+0.16, 6.72+0.21, 12.31+0.25],
[7.69, 6.73, 18.68, 18.26, 7.72, 7.60, 18.73, 18.62, 13.00],
]

Legends = [
'PhysBAM:Compute',
#'PhysBAM:Serialize/Deserialize',
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
    p = plot_bar(RefData[i], ind, bottom, SubColors[i], '', False,
                 i==len(RefData)-1)
    Parts.append(p[0])
bottom = [0] * N
for i in range(0, len(Data)):
    p = plot_bar(Data[i], ind, bottom, Colors[i], '', True,
                 i==len(Data)-1)
    Parts.append(p[0])
pos = sum(i[0] for i in RefData)
# plt.annotate('PhysBAM', xy=(width,pos), xytext=(width,pos+5),
#              fontsize='xx-small',
#              arrowprops=dict(facecolor='black',width=1.2,frac=.3),
#              horizontalalignment='center',
#              verticalalignment='center')
pos = sum(i[0] for i in Data)
# plt.annotate('Nimbus', xy=(0,pos), xytext=(0,pos+5),
#              fontsize='xx-small',
#              arrowprops=dict(facecolor='black',width=1.2,frac=.3),
#              horizontalalignment='center',
#              verticalalignment='center')


ytop = 4
ticks=['worker 1','worker 2','worker 3','worker 4','worker 5', 'worker 6', 'worker 7', 'worker 8', 'Average']
plt.xticks(ind+width/2., ticks)
plt.yticks(np.linspace(0,ytop*10,ytop+1))
# plt.xlabel('Number of workers')
plt.ylabel('Time(s)')

margin = 0.2
plt.xlim(-0.2-margin,N-0.5+.6+margin)
plt.ylim(0,35)
plt.legend(reversed(Parts), reversed(Legends),
           ncol=2, loc=3, bbox_to_anchor=(0.,1.02,1.,.102),
           borderaxespad=0., mode='expand',
           fontsize='small', frameon=False)

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.grid(axis='y')
plt.gca().get_xaxis().set_tick_params(top='off')
plt.tight_layout(pad=0.1, rect=(0,0,1,0.8))

plt.show()
# plt.savefig('../figs/weak_scale.pdf')
