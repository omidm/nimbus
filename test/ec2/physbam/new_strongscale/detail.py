#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

import sys

def mean(x):
    return sum(x)/len(x)

d = sys.argv[1]

Colors = []
Colors.append('#8dd3c7')
Colors.append('#fdc086')
Colors.append('#ffffb3')
Colors.append('#bebada')

a=[]
b=[]
c=[]
rank = 1
while True:
    packtime = 0
    synctime = 0
    try:
        f = open('{}mpi{}.log'.format(d, rank))
    except IOError as e:
        N = rank - 1
        break
    steps = 0
    rank += 1
    for line in f:
        if 'Pack' in line or 'Prepare' in line or 'Unpack' in line:
            packtime += float(line.split()[1])
        elif 'Wait' in line or 'Reduction' in line:
            synctime += float(line.split()[1])
        elif 'Simulation' in line:
            calctime = float(line.split()[1]) - packtime - synctime
        elif 'Substep' in line:
            steps += 1
    plt.bar([rank-1-0.5], calctime/steps, width=1, color=Colors[0])
    plt.bar([rank-1-0.5], packtime/steps,
            bottom=calctime/steps, width=1, color=Colors[1])
    plt.bar([rank-1-0.5], synctime/steps,
            bottom=(packtime+calctime)/steps,
            width=1, color=Colors[3])
    a.append(packtime/steps)
    b.append(synctime/steps)
    c.append(calctime/steps)
plt.xticks(range(8, N+1, 8),range(8, N+1, 8))
plt.xlabel('MPI processor rank')
plt.ylabel('time(s) per iteration')
#plt.ylim(0,50)
plt.xlim(0+.5,N+.5)
plt.title('#Worker {} Substeps {} Frames 3\n'
          'calc {:.2f} pack {:.2f} sync {:.2f}'.format(
              N, steps, mean(c), mean(a), mean(b)))
plt.savefig('w1.pdf')
#plt.show()

