#!/usr/bin/env python 
import copy
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import numpy as np
import pylab

##############################################################################
#                                 PARAMETERS                                 #
##############################################################################
num_threads = 8
num_workers = (int(sys.argv[2]) - int(sys.argv[1]) + 1)
# set denom_ave to number of iterations or frames to get average value instead
# of aggregate
denom_ave = 55

##############################################################################
#                            GROUP INTO CATGORIES                            #
##############################################################################

times_grp = {}

# NOTE: to add categories in plot, edit the categories here.

times_grp['copy']                = {}
times_grp['copy']['cache_total'] = []
times_grp['copy']['block']       = []
times_grp['copy']['lock']        = []
times_grp['copy']['wfc']         = []
times_grp['copy']['mapping']     = []

times_grp['cv']                = {}
times_grp['cv']['cache_total'] = []
times_grp['cv']['block']       = []
times_grp['cv']['lock']        = []
times_grp['cv']['wfc']         = []
times_grp['cv']['rtc']         = []
times_grp['cv']['mapping']     = []

times_grp['cs']                = {}
times_grp['cs']['cache_total'] = []
times_grp['cs']['block']       = []
times_grp['cs']['lock']        = []
times_grp['cs']['wfc']         = []
times_grp['cs']['rtc']         = []
times_grp['cs']['mapping']     = []

# group parsed numbers into plotting categories for each file
for i in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    for k1 in times_grp.keys():
        for k2 in times_grp[k1].keys():
            times_grp[k1][k2].append(0)
    fname = str(i) + "_parse.txt"
    f = open(fname)
    for line in f:
        line = line.strip("\t\n\r\f\v ,}")

        #######################################################################

        # NOTE: to add categories in plot, edit the parser here.

        if "stage" in line and "Copy job" in line:
            times_grp['copy']['cache_total'][-1] += float(line.split()[-1])

        elif "block" in line and "Copy job" in line:
            times_grp['copy']['block'][-1] += float(line.split()[-1])

        elif " lock" in line and  "Copy job" in line:
            times_grp['copy']['lock'][-1] += float(line.split()[-1])

        elif " mapping" in line and  "Copy job" in line:
            times_grp['copy']['mapping'][-1] += float(line.split()[-1])

        elif "wfc" in line and "Copy job" in line:
            times_grp['copy']['wfc'][-1] += float(line.split()[-1])

        #######################################################################

        # NOTE: to add categories in plot, edit the parser here.

        elif "GAV stage" in line or "WIV stage" in line:
            times_grp['cv']['cache_total'][-1] += float(line.split()[-1])
        elif "GAS stage" in line or "WIS stage" in line:
            times_grp['cs']['cache_total'][-1] += float(line.split()[-1])

        elif "GAV block" in line or "WIV block" in line:
            times_grp['cv']['block'][-1] += float(line.split()[-1])
        elif "GAS block" in line or "WIS block" in line:
            times_grp['cs']['block'][-1] += float(line.split()[-1])

        elif "GAV lock" in line or "WIV lock" in line:
            times_grp['cv']['lock'][-1] += float(line.split()[-1])
        elif "GAS lock" in line or "WIS lock" in line:
            times_grp['cs']['lock'][-1] += float(line.split()[-1])

        elif "GAV mapping" in line or "WIV mapping" in line:
            times_grp['cv']['mapping'][-1] += float(line.split()[-1])
        elif "GAS mapping" in line or "WIS mapping" in line:
            times_grp['cs']['mapping'][-1] += float(line.split()[-1])

        elif "GAV wfc" in line or "WIV wfc" in line:
            times_grp['cv']['wfc'][-1] += float(line.split()[-1])
        elif "GAS wfc" in line or "WIS wfc" in line:
            times_grp['cs']['wfc'][-1] += float(line.split()[-1])

        elif "GAV rtc" in line or "WIV rtc" in line:
            times_grp['cv']['rtc'][-1] += float(line.split()[-1])
        elif "GAS rtc" in line or "WIS rtc" in line:
            times_grp['cs']['rtc'][-1] += float(line.split()[-1])

        #######################################################################

# average over all threads
for k1 in times_grp.keys():
    for k2 in times_grp[k1].keys():
        times_grp[k1][k2] = map(lambda x : x/(num_threads*denom_ave), times_grp[k1][k2])

# aggregated, uncategorized into other category
times_grp_calc = {}
times_grp_calc['copy'] = [0] * num_workers
times_grp_calc['cv'] = [0] * num_workers
times_grp_calc['cs'] = [0] * num_workers
for k1 in times_grp.keys():
    for k2 in times_grp[k1].keys():
        if k2 != 'cache_total':
            for i in range(0, num_workers):
                times_grp_calc[k1][i] += times_grp[k1][k2][i]
for k in times_grp.keys():
    times_grp[k]['other'] = []
    for i in range(0, num_workers):
        times_grp[k]['other'].append(times_grp[k]['cache_total'][i] - \
                                     times_grp_calc[k][i])

##############################################################################
#                     PLOT CATGEORIES INTO STACKED PLOTS                     #
##############################################################################

# NOTE: to add categories in plot, edit number of categories here.
copy_num = 5
cv_num = 6
cs_num = 6
num_categories = copy_num + cv_num + cs_num

# plot properties
ind = np.arange(num_workers)
width = 4 * 1 / (num_categories + 1.0)
pylab.rcParams['xtick.major.pad']='20'
pylab.rcParams['ytick.major.pad']='20'
font = {'family' : 'sans-serif',
        'size'   : 12,
        'weight' : 'bold'}
matplotlib.rc('font', **font)
fig, ax = plt.subplots()
fig.patch.set_facecolor('white')

# plot data
def addTimes(a, b):
    times = a
    for idx, val in enumerate(b):
        times[idx] += val
    return times

groups = []
for i in range(1, num_workers+1):
    groups.append("W " + str(i))

# NOTE: to add categories in plot, edit number of *num, categories, dataList here.
categories = [ 'Copy Block', 'Copy Lock', 'Copy WriteFromCache', \
               'Copy Mapping', 'Copy Other', \
               'Comp Var Block', 'Comp Var Lock', 'Comp Var WriteFromCache', \
               'Comp Var ReadToCache', 'Comp Var Mapping', 'Comp Var Other', \
               'Comp Struct Block', 'Comp Struct Lock', 'Comp Struct WriteFromCache', \
               'Comp Struct ReadToCache', 'Comp Struct Mapping', 'Comp Struct Other' ]
dataList = [ times_grp['copy']['block'], times_grp['copy']['lock'], times_grp['copy']['wfc'], \
             times_grp['copy']['mapping'], times_grp['copy']['other'], \
             times_grp['cv']['block'], times_grp['cv']['lock'], \
             times_grp['cv']['wfc'], times_grp['cv']['rtc'], \
             times_grp['cv']['mapping'], times_grp['cv']['other'], \
             times_grp['cs']['block'], times_grp['cs']['lock'], \
             times_grp['cs']['wfc'], times_grp['cs']['rtc'], \
             times_grp['cs']['mapping'], times_grp['cs']['other'] ]

times = {}
for i in range(len(categories)):
    times[categories[i]] = dataList[i]
timesCopy = copy.deepcopy(times)
bottoms = {}
bottoms[categories[0]] = [0] * num_workers
for i in range(1, num_categories):
    bottoms[categories[i]] = addTimes(timesCopy[categories[i-1]], \
                                      bottoms[categories[i-1]])

# colors and hashing
colors = []
colors.append('#8dd3c7')
colors.append('#ffffb3')
colors.append('#bebada')
colors.append('#fb8072')
colors_use = [colors[0]]*copy_num + \
             [colors[2]]*cv_num + [colors[3]]*cs_num

# NOTE: to add categories in plot, edit hatch patterns here.
hatch_block = "x"
hatch_lock = "*"
hatch_write = "/"
hatch_read  = "\\"
hatch_map = "|"
hatch_use  = [ hatch_block, hatch_lock, hatch_write, hatch_map, "", \
               hatch_block, hatch_lock, hatch_write, hatch_read, hatch_map, "", \
               hatch_block, hatch_lock, hatch_write, hatch_read, hatch_map, "" ]

# bars
rects = {}
for i in range(0, num_categories):
    rects[categories[i]] = ax.bar(ind, times[categories[i]], 1.5*width, \
                                  bottom=bottoms[categories[i]], \
                                  color=colors_use[i], hatch=hatch_use[i]) 

# label bars with numbers (values) for selected categories
def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text((rect.get_x() + rect.get_width() * 1.05), \
                 rect.get_y() + 0.5*height - 0.05, '%0.2f'%(height), \
                 ha='left', va='bottom', fontsize='8')
#label_list = { 'LC WriteFromCache', 'LC Block', \
#               'RC WriteFromCache', 'RC Block', \
#               'Comp Var WriteFromCache', 'Comp Var ReadToCache', 'Comp Var Block' }
#for i in label_list:
#    autolabel(rects[i])

# label and ticks
ax.set_ylabel('Time (s)')
ax.set_title('Cache Times - Scale 256, 64 Partitions (Uniform), 3 Frames')
ax.set_xticks(ind + 0.75*width)
ax.set_xticklabels( groups )
ax.set_yticks(np.arange(0, 5.1, 0.5))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.yaxis.label.set_weight('bold')
ax.yaxis.labelpad = 6

# legend
# NOTE: to add categories in plot, edit legend here.
blockBars = ax.bar(0, 0, 0, 0, color='#FFFFFF', hatch=hatch_block)
lockBars = ax.bar(0, 0, 0, 0, color='#FFFFFF', hatch=hatch_lock)
writeBars = ax.bar(0, 0, 0, 0, color='#FFFFFF', hatch=hatch_write)
readBars = ax.bar(0, 0, 0, 0, color='#FFFFFF', hatch=hatch_read)
mapBars = ax.bar(0, 0, 0, 0, color='#FFFFFF', hatch=hatch_map)
leg = ax.legend( (rects['Comp Struct Other'][0], rects['Comp Var Other'][0], \
                  rects['Copy Other'][0], \
                  blockBars[0], lockBars[0], readBars[0], writeBars[0], mapBars[0]), \
                 ('Compute Job - AppStruct', 'Compute Job - AppVar', 'Copy Jobs', \
                  'Block Time', 'Lock Time', \
                  'Read to Cache', 'Write from Cache', 'Edit Mapping'), \
                 labelspacing = 0.0, borderpad =  0.2, loc = 9,
                 ncol = 2)
leg.draggable()
for label in leg.get_texts():
    label.set_fontsize(8)

# scale
x1,x2,y1,y2 = plt.axis()
plt.axis((x1-0.5,x2+0.5,0,5.5))

# grid
plt.grid()

# save figure
plt.savefig("cache_times.pdf")
