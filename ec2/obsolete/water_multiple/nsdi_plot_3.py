#!/usr/bin/python
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import rc

# Configure matplot to understand tex grammer.
# rc('text', usetex=True)
rc('font', size=24)

N = 12
ind = np.arange(N) + 1
width = 0.5
shift = 0.2

Legends = [
'MPI',
'Nimbus /w  CT + WT'
]

n_colors = 4
color_map = cmx.ScalarMappable(
        colors.Normalize(vmin=0, vmax=n_colors+1), cmap='Greens')
Colors = [color_map.to_rgba(i) for i in range(n_colors, -1, -1)]

color_map = cmx.ScalarMappable(
        colors.Normalize(vmin=0, vmax=n_colors+1), cmap='Reds')
RColors = [color_map.to_rgba(i) for i in range(n_colors, -1, -1)]


nimbus = [190, 99.62, 146, 54.32, 38.41, 37.31, 36.99, 37, 37.3, 36.68, 36.66, 36.70]
mpi    = [31, 30.4, 31.1, 30.68, 31.2, 31, 31.1, 31.4, 31.4, 31.4, 31.1, 31.5]


Parts = []

line = plt.bar(ind + shift, nimbus, width, color=Colors[3])
Parts.append(line[0])

line = plt.bar(ind, mpi, width, color=RColors[0])
Parts.append(line[0])



# ticks=['(8, 64, 512$^3$)','(64, 512, 1024$^3$)']
plt.xticks(ind)
plt.xlim([1, 12])
plt.xlabel('Iteration number', fontsize='medium')
plt.ylabel('Iteration length (seconds)', fontsize='medium')
# plt.ylim([0, 210])
# # plt.set_yscale('log')

plt.legend(reversed(Parts), Legends,
           ncol=1, loc='best',
           # mode='expand',
           fontsize='medium', frameon=False)

plt.show()
