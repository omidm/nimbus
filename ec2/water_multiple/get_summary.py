#!/usr/bin/python

# ./get_summary.py [directory_name] [number of workers]

import sys

def parse_line(line):
    items = line.split()
    return items[3], float(items[5])

s1=s2=s3=s8=s5=0
stats = list()
N = int(sys.argv[2])
for rank in range(1, N+1):
    f = open('{}/{}_time_per_thread.txt'.format(sys.argv[1],rank), 'r')

    stat = dict()

    for line in f:
        key, time = parse_line(line)
        if key in stat:
            stat[key] += time
        else:
            stat[key] = time
    stats.append(stat)

    steps=stat['kCalculateDt'] /8.
    s1+=steps
    comp=(stat['kExecuteComputationJob']-stat['kAssemblingCache']) /steps/8.
    s2+=comp
    cache=(stat['kExecuteCopyJob']+stat['kAssemblingCache']) /steps/8.
    s3+=cache
    block=stat['kSumCyclesRun']/steps/8-comp-cache+stat['kSumCyclesBlock'] /8/steps
    s8+=block
    idle=(stat['kSumCyclesTotal']-stat['kSumCyclesBlock']-stat['kSumCyclesRun']) /steps/8
    s5+=idle
    print '{:8.0f}: {:8.0f} {:8.2f} {:8.2f} {:8.2f} {:8.2f}'.format(
            rank, steps, comp, cache, block, idle)
print '   steps  compute cache/copy  block    idle:'
print '{:8.0f} {:8.2f} {:8.2f} {:8.2f} {:8.2f}'.format(s1/N, s2/N, s3/N, s8/N, s5/N)
