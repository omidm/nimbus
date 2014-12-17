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
num_categories = 14
# set denom_ave to number of iterations or frames to get average value instead
# of aggregate
denom_ave = 55

##############################################################################
#                            GROUP INTO CATGORIES                            #
##############################################################################

lc_cache_total = []
lc_block       = []
lc_wfc         = []

rc_cache_total = []
rc_block       = []
rc_wfc         = []

cv_cache_total = []
cv_block       = []
cv_wfc         = []
cv_rtc         = []

cs_cache_total = []
cs_block       = []
cs_wfc         = []
cs_rtc         = []

# group parsed numbers into plotting categories for each file
for i in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    lc_cache_total.append(0)
    lc_block.append(0)
    lc_wfc.append(0)
    rc_cache_total.append(0)
    rc_block.append(0)
    rc_wfc.append(0)
    cv_cache_total.append(0)
    cv_block.append(0)
    cv_wfc.append(0)
    cv_rtc.append(0)
    cs_cache_total.append(0)
    cs_block.append(0)
    cs_wfc.append(0)
    cs_rtc.append(0)
    fname = str(i) + "_parse.txt"
    f = open(fname)
    for line in f:
        line = line.strip("\t\n\r\f\v ,}")
        if "stage LC job" in line:
            lc_cache_total[-1] += float(line.split()[-1])
        elif ("stage RCS job" in line or "stage RCR job" in line):
            rc_cache_total[-1] += float(line.split()[-1])
        elif "block LC job" in line:
            lc_block[-1] += float(line.split()[-1])
        elif ("block RCS job" in line or "block RCR job" in line):
            rc_block[-1] += float(line.split()[-1])
        elif " lock LC job" in line:
            lc_block[-1] += float(line.split()[-1])
        elif (" lock RCS job" in line or " lock RCR job" in line):
            rc_block[-1] += float(line.split()[-1])
        elif "pdata" in line and "LC job" in line:
            lc_wfc[-1] += float(line.split()[-1])
        elif "pdata" in line and ("RCR job" in line or "RCS job" in line):
            rc_wfc[-1] += float(line.split()[-1])
        elif "GAV stage" in line or "WIV stage" in line:
            cv_cache_total[-1] += float(line.split()[-1])
        elif "GAS stage" in line or "WIS stage" in line:
            cs_cache_total[-1] += float(line.split()[-1])
        elif "GAV block" in line or "WIV block" in line:
            cv_block[-1] += float(line.split()[-1])
        elif "GAS block" in line or "WIS block" in line:
            cs_block[-1] += float(line.split()[-1])
        elif "GAV lock" in line or "WIV lock" in line:
            cv_block[-1] += float(line.split()[-1])
        elif "GAS lock" in line or "WIS lock" in line:
            cs_block[-1] += float(line.split()[-1])
        elif "GAV wfc" in line or "WIV wfc" in line:
            cv_wfc[-1] += float(line.split()[-1])
        elif "GAS wfc" in line or "WIS wfc" in line:
            cs_wfc[-1] += float(line.split()[-1])
        elif "GAV rtc" in line or "WIV rtc" in line:
            cv_rtc[-1] += float(line.split()[-1])
        elif "GAS rtc" in line or "WIS rtc" in line:
            cs_rtc[-1] += float(line.split()[-1])

# average over all threads
lc_cache_total = map(lambda x : x/(num_threads*denom_ave), lc_cache_total)
lc_block = map(lambda x : x/(num_threads*denom_ave), lc_block)
lc_wfc = map(lambda x : x/(num_threads*denom_ave), lc_wfc)
rc_cache_total = map(lambda x : x/(num_threads*denom_ave), rc_cache_total)
rc_block = map(lambda x : x/(num_threads*denom_ave), rc_block)
rc_wfc = map(lambda x : x/(num_threads*denom_ave), rc_wfc)
cv_cache_total = map(lambda x : x/(num_threads*denom_ave), cv_cache_total)
cv_block = map(lambda x : x/(num_threads*denom_ave), cv_block)
cv_wfc = map(lambda x : x/(num_threads*denom_ave), cv_wfc)
cv_rtc = map(lambda x : x/(num_threads*denom_ave), cv_rtc)
cs_cache_total = map(lambda x : x/(num_threads*denom_ave), cs_cache_total)
cs_block = map(lambda x : x/(num_threads*denom_ave), cs_block)
cs_wfc = map(lambda x : x/(num_threads*denom_ave), cs_wfc)
cs_rtc = map(lambda x : x/(num_threads*denom_ave), cs_rtc)

# aggregated, uncategorized into other category
lc_other = []
rc_other = []
cv_other = []
cs_other = []
for i in range(0, num_workers):
    lc_other.append(lc_cache_total[i] - lc_block[i] - lc_wfc[i])
    rc_other.append(rc_cache_total[i] - rc_block[i] - rc_wfc[i])
    cv_other.append(cv_cache_total[i] - cv_block[i] - cv_wfc[i] - cv_rtc[i])
    cs_other.append(cs_cache_total[i] - cs_block[i] - cs_wfc[i] - cs_rtc[i])

##############################################################################
#                     PLOT CATGEORIES INTO STACKED PLOTS                     #
##############################################################################

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
categories = [ 'LC Block', 'LC WriteFromCache', 'LC Other', \
               'RC Block', 'RC WriteFromCache', 'RC Other', \
               'Comp Var Block', 'Comp Var WriteFromCache', \
               'Comp Var ReadToCache', 'Comp Var Other', \
               'Comp Struct Block', 'Comp Struct WriteFromCache', \
               'Comp Struct ReadToCache', 'Comp Struct Other' ]
dataList = [ lc_block, lc_wfc, lc_other, \
             rc_block, rc_wfc, rc_other, \
             cv_block, cv_wfc, cv_rtc, cv_other, \
             cs_block, cs_wfc, cs_rtc, cs_other ]
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
colors_use = [colors[0]]*3 + [colors[1]]*3 + [colors[2]]*4 + [colors[3]]*4
hatch_block = "x"
hatch_write = "/"
hatch_read  = "\\"
hatch_use  = [ hatch_block, hatch_write, "", \
               hatch_block, hatch_write, "", \
               hatch_block, hatch_write, hatch_read, "", \
               hatch_block, hatch_write, hatch_read, "" ]

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
label_list = { 'LC WriteFromCache', 'LC Block', \
               'RC WriteFromCache', 'RC Block', \
               'Comp Var WriteFromCache', 'Comp Var ReadToCache', 'Comp Var Block' }
for i in label_list:
    autolabel(rects[i])

# label and ticks
ax.set_ylabel('Time (s)')
ax.set_title('Cache Times - Scale 256, 64 Partitions (Uniform), 3 Frames')
ax.set_xticks(ind + 0.75*width)
ax.set_xticklabels( groups )
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.yaxis.label.set_weight('bold')
ax.yaxis.labelpad = 6

# legend
blockBars = ax.bar(0, 0, 0, 0, color='#FFFFFF', hatch=hatch_block)
writeBars = ax.bar(0, 0, 0, 0, color='#FFFFFF', hatch=hatch_write)
readBars = ax.bar(0, 0, 0, 0, color='#FFFFFF', hatch=hatch_read)
leg = ax.legend( (rects['Comp Struct Other'][0], rects['Comp Var Other'][0], \
                  rects['RC Other'][0], rects['LC Other'][0], \
                  blockBars[0], readBars[0], writeBars[0]), \
                 ('Compute Job - AppStruct', 'Compute Job - AppVar', \
                  'Remote Copy', 'Local Copy', 'Blocked Time', \
                  'Read to Cache', 'Write from Cache'), \
                 labelspacing = 0.0, borderpad =  0.2, loc = 9,
                 ncol = 2)
leg.draggable()
for label in leg.get_texts():
    label.set_fontsize(8)

# scale
x1,x2,y1,y2 = plt.axis()
plt.axis((x1-0.5,x2+0.5,0,5.5))

# save figure
plt.savefig("cache_times.pdf")
