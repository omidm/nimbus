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
ind = np.arange(N - 1)
ind = np.append(ind, [N - 0.5])
width = 0.4


### Data, btaching size 0f 10 only for 1024 but still has bottleneck job 
Data = [
# run time
[
[11.90,12.20], # avg(run)
],

# synch time
[
[12.09, 19.42], # avg(block) + [avg(idle) - min_worker(idle)]
],

# idle time
[
# [ 0.67, 3.38], # controller(AC+DMQ+TI)
# [ 0.14, 0.83], # parent_exec
# [ 0.52, 0.71]  # min_worker(idle) - controller(AC+DMQ+TI) - parent_exec
[1.33, 4.92]
]

]


RefData = [
# run time
[
[11.85, 12.42], # avg(run)
],

# synch time
[
[12.88, 19.27], # avg(block)
]

]



Legends = [
'MPI:Compute',
'MPI:Comm',
'Nimbus:Compute',
'Nimbus:Comm',
'Nimbus:Scheduler',
]

n_colors = 3
color_map = cmx.ScalarMappable(
        colors.Normalize(vmin=0, vmax=n_colors+1), cmap='Greens')
Colors = [color_map.to_rgba(i) for i in range(n_colors, -1, -1)]

color_map = cmx.ScalarMappable(
        colors.Normalize(vmin=0, vmax=n_colors+1), cmap='Reds')
SubColors = [color_map.to_rgba(i) for i in range(n_colors, -1, -1)]

Hatch = ['', '\\\\', '//']





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
      if (num > 1.5):
        plt.text(rect.get_x() + rect.get_width()/2.,
                 rect.get_y()+ rect.get_height()/2.,
                 '{:.1f}'.format(num),
                 ha='center', va='center',
                 fontsize='small')
      if add_sum:
          plt.text(rect.get_x() + rect.get_width()/2,
                   rect.get_y()+ rect.get_height() + 1,
                   '{:.1f}'.format(num+bottom[i]),
                   ha='center', va='center',
                   fontsize='small')
          i += 1
    for i in range(0, len(bottom)):
        bottom[i] += bar_data[i]
    return p

Parts = []
bottom = [0] * N
for i in range(0, len(RefData)):
    for j in range(0, len(RefData[i])):
      p = plot_bar(RefData[i][j], ind, bottom, SubColors[i], Hatch[j], True, (j==len(RefData[i])-1) and (i==len(RefData)-1))
      Parts.append(p[0])


bottom = [0] * N
for i in range(0, len(Data)):
    for j in range(0, len(Data[i])):
      p = plot_bar(Data[i][j], ind, bottom, Colors[i], Hatch[j], False, (j==len(Data[i])-1) and (i==len(Data)-1))
      Parts.append(p[0])







ticks=['(8, 64, 512$^3$)','(64, 512, 1024$^3$)']
plt.xticks(ind+width/2., ticks, fontsize='medium')
plt.xlabel('(#workers, #cores, #cells)', fontsize='large')
plt.yticks(np.linspace(0,60,7))
plt.ylabel('Time (seconds)', fontsize='large')

margin = 0.2
plt.xlim(-0.2-margin,N-0.5+.6+margin)
# plt.ylim(0,35)
fl = plt.legend(reversed(Parts[0:2]), reversed(Legends[0:2]),
           ncol=1, loc=2,
           # bbox_to_anchor=(0.,1,1,0),
           borderaxespad=0.,
           # mode='expand',
           fontsize='large', frameon=False)

plt.gca().add_artist(fl)


plt.legend(reversed(Parts[2:5]), reversed(Legends[2:5]),
           ncol=1, loc=1,
           # bbox_to_anchor=(0.,1,1,0),
           borderaxespad=0.,
           # mode='expand',
           fontsize='large', frameon=False)




# plt.gca().spines['top'].set_visible(False)
# plt.gca().spines['right'].set_visible(False)
# plt.grid(axis='y')
# plt.gca().get_xaxis().set_tick_params(top='off')
# plt.tight_layout(pad=0.1, rect=(0,0,1,0.8))

plt.show()
# plt.savefig('../figs/weak_scale.pdf')




