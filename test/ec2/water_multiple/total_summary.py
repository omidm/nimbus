#!/usr/bin/python

# ./get_summary.py [directory_name] [number of workers]

import sys
import argparse

def parse_line(line):
    items = line.split()
    return items[3], float(items[5])



## Parse the command line arguments ##
parser = argparse.ArgumentParser(description='Process log files.')
parser.add_argument(
    "-d", "--dir_path",
    dest="dirpath",
    help="input dorectory path")
parser.add_argument(
    "-wn", "--worker_num",
    dest="workernum",
    default=8,
    help="worker num")
parser.add_argument(
    "-cn", "--core_num",
    dest="corenum",
    default=8,
    help="core num per worker")

args = parser.parse_args()

d  = args.dirpath 
N  = int(args.workernum)
CN = int(args.corenum)







print '--------------------------------------------------------------------------------------'
print '   worker    steps  compute   c/copy      run      block     idle    synch     total'
print '--------------------------------------------------------------------------------------'

s1=s2=s3=s8=s5=0
stats = list()

for rank in range(1, N+1):
    f = open('{}/{}_time_per_thread.txt'.format(d, rank), 'r')

    stat = dict()

    for line in f:
        key, time = parse_line(line)
        if key in stat:
            stat[key] += time
        else:
            stat[key] = time
    stats.append(stat)

    steps=stat['kCalculateDt'] /CN
    s1+=steps
    comp=(stat['kExecuteComputationJob']-stat['kAssemblingCache']) /steps/CN
    s2+=comp
    cache=(stat['kExecuteCopyJob']+stat['kAssemblingCache']) /steps/CN
    s3+=cache
    block=stat['kSumCyclesRun']/steps/CN-comp-cache+stat['kSumCyclesBlock'] /CN/steps
    s8+=block
    idle=(stat['kSumCyclesTotal']-stat['kSumCyclesBlock']-stat['kSumCyclesRun']) /steps/CN
    s5+=idle
    print '{:8.0f}: {:8.0f} {:8.2f} {:8.2f} {:8.2f}   {:8.2f} {:8.2f} {:8.2f}  {:8.2f}'.format(
            rank, steps, comp, cache, comp+cache, block, idle, block+idle, comp+cache+block+idle)

print '--------------------------------------------------------------------------------------'
print ' Average: {:8.0f} {:8.2f} {:8.2f} {:8.2f}   {:8.2f} {:8.2f} {:8.2f}  {:8.2f}'.format(s1/N, s2/N, s3/N, (s2+s3)/N, s8/N, s5/N, (s8+s5)/N ,(s2+s3+s8+s5)/N)
print '--------------------------------------------------------------------------------------'
