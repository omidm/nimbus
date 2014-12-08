    #!/usr/bin/env python
import sys
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib
import pylab
import copy
from matplotlib import colors

def addTimes(a, b):
    times = a
    for idx, val in enumerate(b):
        times[idx] += val
    return times

categories = ['SyncData Check LC', 'SyncData PullData LC', 'SyncData Other LC', 'InvalidateMappings Check LC', 'InvalidateMappings Other LC', \
    'SyncData Check RC', 'SyncData PullData RC', 'SyncData Other RC', 'InvalidateMappings Check RC', 'InvalidateMappings Other RC', \
    'GetAppStruct Check', 'GetAppStruct ReadFromCache', 'GetAppStruct WriteToCache', 'GetAppStruct Other', \
    'GetAppVar Check', 'GetAppVar ReadFromCache', 'GetAppVar WriteToCache', 'GetAppVar Other']

syncDataCheckLC = [3.4370462894439697, 1.1585783958435059, 2.0288302898406982, 2.022611141204834, 2.250786781311035, 3.2912967205047607, 4.543894290924072, 2.4757838249206543]
syncDataTotalLC = [57.8506383895874, 54.58799147605896, 54.72220993041992, 57.459338903427124, 55.65745162963867, 60.34208846092224, 59.843356132507324, 57.237120389938354]
syncDataPullDataLC = [45.469523668289185, 44.58235025405884, 43.61118841171265, 46.2102484703064, 43.967732667922974, 47.79056739807129, 46.276620864868164, 45.763588190078735]
syncDataOtherLC = []
for i in range(8):
    syncDataOtherLC.append(syncDataTotalLC[i] - syncDataCheckLC[i] - syncDataPullDataLC[i])

syncDataCheckRC = [44.44410300254822, 17.266765356063843, 18.802417993545532, 18.17496371269226, 16.663784503936768, 28.70637559890747, 29.766319274902344, 47.451478719711304]
syncDataTotalRC = [105.53758668899536, 76.31077122688293, 80.42673110961914, 72.3644449710846, 73.46345901489258, 87.78239798545837, 87.86790823936462, 110.6665210723877]
syncDataPullDataRC = [52.86698794364929, 51.200355052948, 54.037930488586426, 46.26099252700806, 49.108949422836304, 50.49369525909424, 49.487077951431274, 54.783613443374634]
syncDataOtherRC = []
for i in range(8):
    syncDataOtherRC.append(syncDataTotalRC[i] - syncDataCheckRC[i] - syncDataPullDataRC[i])

invalidateMappingsCheckLC = [0.3688807487487793, 0.3682892322540283, 0.38413214683532715, 0.386946439743042, 0.3636164665222168, 0.38910770416259766, 0.3682279586791992, 0.37645483016967773]
invalidateMappingsTotalLC = [2.6446011066436768, 2.620304822921753, 2.7446165084838867, 2.8067164421081543, 2.6702823638916016, 2.7701499462127686, 2.673391580581665, 2.7565059661865234]
invalidateMappingsOtherLC = []
for i in range(8):
    invalidateMappingsOtherLC.append(invalidateMappingsTotalLC[i] - invalidateMappingsCheckLC[i])

invalidateMappingsCheckRC = [0.4514350891113281, 0.43439197540283203, 0.45607995986938477, 0.7300543785095215, 0.4472770690917969, 0.45235300064086914, 0.4342002868652344, 0.4569213390350342]
invalidateMappingsTotalRC = [4.38065242767334, 4.265584945678711, 4.506795644760132, 6.771206855773926, 4.205585479736328, 4.619951009750366, 4.323426246643066, 4.546139478683472]
invalidateMappingsOtherRC = []
for i in range(8):
    invalidateMappingsOtherRC.append(invalidateMappingsTotalRC[i] - invalidateMappingsCheckRC[i])

getAppStructCheck = [0.10576105117797852, 0.10277223587036133, 0.10534930229187012, 0.10796451568603516, 0.10245394706726074, 0.11881208419799805, 0.10532450675964355, 0.10748672485351562]
getAppStructReadFromCache = [8.745710849761963, 6.30970573425293, 7.956263065338135, 18.56783390045166, 5.725188732147217, 6.25141978263855, 6.734279155731201, 8.26154899597168]
getAppStructWriteToCache = [2.1199965476989746, 1.5326061248779297, 1.7818207740783691, 2.711050033569336, 1.4223268032073975, 1.7194252014160156, 1.8287007808685303, 1.9950296878814697]
getAppStructTotal = [11.717333316802979, 8.706499099731445, 10.552387237548828, 22.748674154281616, 7.926727056503296, 8.96130633354187, 9.518373966217041, 11.091608047485352]
getAppStructOther = []
for i in range(8):
    getAppStructOther.append(getAppStructTotal[i] - getAppStructCheck[i] - getAppStructWriteToCache[i] - getAppStructReadFromCache[i])

getAppVarCheck = [43.307443618774414, 29.886924505233765, 53.57873201370239, 63.49052453041077, 29.613428592681885, 32.32182312011719, 27.691977500915527, 40.92560911178589]
getAppVarReadFromCache = [40.45420861244202, 38.7151825428009, 39.687335729599, 54.8624222278595, 38.32898998260498, 44.243202209472656, 40.319087982177734, 41.04272484779358]
getAppVarWriteToCache = [32.971394777297974, 21.534605503082275, 30.885990858078003, 36.83485817909241, 21.584059953689575, 24.907496452331543, 23.71062660217285, 32.69905877113342]
getAppVarTotal = [134.9231550693512, 107.652348279953, 144.049964427948, 172.86738276481628, 106.97693848609924, 119.58382821083069, 109.15449547767639, 133.5821831226349]
getAppVarOther = []
for i in range(8):
    getAppVarOther.append(getAppVarTotal[i] - getAppVarCheck[i] - getAppVarWriteToCache[i] - getAppVarReadFromCache[i])

groups = ['Worker 1', 'Worker 2', 'Worker 3', 'Worker 4', 'Worker 5', 'Worker 6', 'Worker 7', 'Worker 8']

times = {}

dataList = [syncDataCheckLC, syncDataPullDataLC, syncDataOtherLC, \
    invalidateMappingsCheckLC, invalidateMappingsOtherLC, \
    syncDataCheckRC, syncDataPullDataRC, syncDataOtherRC, \
    invalidateMappingsCheckRC, invalidateMappingsOtherRC, \
    getAppStructCheck, getAppStructReadFromCache, getAppStructWriteToCache, getAppStructOther, \
    getAppVarCheck, getAppVarReadFromCache, getAppVarWriteToCache, getAppVarOther]

for i in range(len(categories)):
    times[categories[i]] = dataList[i]

print times

num_workers = 8
num_categories = 12

ind = np.arange(num_workers)
width = 4 * 1 / (num_categories + 1.0)

pylab.rcParams['xtick.major.pad']='20'
pylab.rcParams['ytick.major.pad']='20'

font = {'family' : 'sans-serif',
        'size'   : 20,
        'weight' : 'bold'}

matplotlib.rc('font', **font)


fig, ax = plt.subplots()
fig.patch.set_facecolor('white')

timesCopy = copy.deepcopy(times)

bottom2 = timesCopy[categories[0]]
print bottom2 
bottom3 = addTimes(timesCopy[categories[1]], bottom2)
print bottom3
bottom4 = addTimes(timesCopy[categories[2]], bottom3)
print bottom4
bottom5 = addTimes(timesCopy[categories[3]], bottom4)
print bottom5
bottom6 = addTimes(timesCopy[categories[4]], bottom5)
print bottom6
bottom7 = addTimes(timesCopy[categories[5]], bottom6)
print bottom7
bottom8 = addTimes(timesCopy[categories[6]], bottom7)
print bottom8
bottom9 = addTimes(timesCopy[categories[7]], bottom8)
print bottom9
bottom10 = addTimes(timesCopy[categories[8]], bottom9)
print bottom10
bottom11 = addTimes(timesCopy[categories[9]], bottom10)
print bottom11
bottom12 = addTimes(timesCopy[categories[10]], bottom11)
print bottom12
bottom13 = addTimes(timesCopy[categories[11]], bottom12)
print bottom13
bottom14 = addTimes(timesCopy[categories[12]], bottom13)
print bottom14
bottom15 = addTimes(timesCopy[categories[13]], bottom14)
print bottom15
bottom16 = addTimes(timesCopy[categories[14]], bottom15)
print bottom16
bottom17 = addTimes(timesCopy[categories[15]], bottom16)
print bottom17
bottom18 = addTimes(timesCopy[categories[16]], bottom17)
print bottom18

colors = []
colors.append('#8dd3c7')
# colors.append('#80b1d3')
colors.append('#ffffb3')
colors.append('#bebada')
colors.append('#fb8072')

rects1 = ax.bar(ind + 3*width, times[categories[0]], width, color=colors[0], hatch="//") 
rects2 = ax.bar(ind + 3*width, times[categories[1]], width, bottom=bottom2, color=colors[0], hatch='x') 
rects3 = ax.bar(ind + 3*width, times[categories[2]], width, bottom=bottom3, color=colors[0])
rects4 = ax.bar(ind + 3*width, times[categories[3]], width, bottom=bottom4, color=colors[0], hatch="//")
rects5 = ax.bar(ind + 3*width, times[categories[4]], width, bottom=bottom5, color=colors[0])
rects6 = ax.bar(ind + 3*width, times[categories[5]], width, bottom=bottom6, color=colors[1], hatch='//') 
rects7 = ax.bar(ind + 3*width, times[categories[6]], width, bottom=bottom7, color=colors[1], hatch='x')
rects8 = ax.bar(ind + 3*width, times[categories[7]], width, bottom=bottom8, color=colors[1])
rects9 = ax.bar(ind + 3*width, times[categories[8]], width, bottom=bottom9, color=colors[1], hatch="//")
rects10 = ax.bar(ind + 3*width, times[categories[9]], width, bottom=bottom10, color=colors[1]) 
rects11 = ax.bar(ind + 3*width, times[categories[10]], width, bottom=bottom11, color=colors[2], hatch="//")
rects12 = ax.bar(ind + 3*width, times[categories[11]], width, bottom=bottom12, color=colors[2], hatch="\\")
rects13 = ax.bar(ind + 3*width, times[categories[12]], width, bottom=bottom13, color=colors[2], hatch="x")
rects14 = ax.bar(ind + 3*width, times[categories[13]], width, bottom=bottom14, color=colors[2])
rects15 = ax.bar(ind + 3*width, times[categories[14]], width, bottom=bottom15, color=colors[3], hatch="//")
rects16 = ax.bar(ind + 3*width, times[categories[15]], width, bottom=bottom16, color=colors[3], hatch="\\")
rects17 = ax.bar(ind + 3*width, times[categories[16]], width, bottom=bottom17, color=colors[3], hatch="x")
rects18 = ax.bar(ind + 3*width, times[categories[17]], width, bottom=bottom18, color=colors[3])

checkBars = ax.bar(0, 0, 0, 0, color='#FFFFFF', hatch="//")
readBars = ax.bar(0, 0, 0, 0, color='#FFFFFF', hatch="\\")
writeBars = ax.bar(0, 0, 0, 0, color='#FFFFFF', hatch="x")

ax.set_ylabel('Time (s)')
ax.set_title('Cache Manager Profiles - Scale 256, 64 Partitions (Uniform), 1 Frame')
ax.set_xticks(ind+3.5*width)
ax.set_xticklabels( ('Worker 1', 'Worker 2', 'Worker 3', 'Worker 4', 'Worker 5', 'Worker 6', 'Worker 7', 'Worker 8') )
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

#ax.title.set_fontsize(18)
#ax.yaxis.label.set_fontsize(50)
ax.yaxis.label.set_weight('bold')
ax.yaxis.labelpad = 15

#leg = ax.legend( (rects1[0], rects2[0], rects3[0], rects4[0], rects5[0], rects6[0]), ('Worker Overhead', 'Copy', 'Load', 'Compute', 'Save', 'Idle'), labelspacing = 0.1, loc = 1)

leg = ax.legend( (rects18[0], rects14[0], rects10[0], rects5[0], checkBars[0], readBars[0], writeBars[0]), \
    ('Compute Job - AppVar Request', 'Compute Job - AppStruct Request', 'Remote Copy', 'Local Copy', 'Blocked Time', 'Read to Cache', 'Write from Cache'), \
    labelspacing = 0.0, borderpad =  0.2, loc = 2)


leg.draggable()
for label in leg.get_texts():
    label.set_fontsize(10)

x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,540))

fig = plt.gcf()
#fig.set_size_inches(20.5,11.5)
#fig.savefig('profiles.png',dpi=100)

def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text((rect.get_x()+rect.get_width())+.05*rect.get_width(), rect.get_y() + height - 5, '%0.2f'%(height + rect.get_y()),
                ha='left', va='bottom', fontsize='10')

autolabel(rects1)
autolabel(rects2)
autolabel(rects6)
autolabel(rects7)

autolabel(rects15)
autolabel(rects16)
autolabel(rects17)

plt.show()

