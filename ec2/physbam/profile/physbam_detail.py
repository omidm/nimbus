#!/usr/bin/python

import numpy as np
import sys
import argparse

def mean(x):
    return sum(x)/len(x)

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
parser.add_argument(
    "-i", "--ignore",
    dest="ignore",
    action="store_true",
    help="ignore the write and initialization in results")

args = parser.parse_args()

d  = args.dirpath 
N  = int(args.workernum)
CN = int(args.corenum)

a=[]
b=[]
c=[]

print '--------------------------------------------------------------------------------------'
print '   worker    steps  compute     pack      run                        synch     total'
print '--------------------------------------------------------------------------------------'

for rank in range(1, CN*N + 1):
    packtime   = 0
    synctime   = 0
    calcsubstep    = 0
    init_phase = True
    f = open('{}mpi{}.log'.format(d, rank))

    steps = 0
    for line in f:
        if init_phase and args.ignore:
          if 'Initialize' in line:
            init_phase = False
          continue

        if 'Pack' in line or 'Prepare' in line or 'Unpack' in line:
            packtime += float(line.split()[1])
        elif 'Wait' in line or 'Reduction' in line:
            synctime += float(line.split()[1])
        elif 'CalculateSubstep' in line:
            steps += 1
            calcsubstep += float(line.split()[1])
        elif 'Simulation' in line:
            if args.ignore:
              calctime = calcsubstep - packtime - synctime
            else:
              calctime = float(line.split()[1]) - packtime - synctime

    a.append(packtime/steps)
    b.append(synctime/steps)
    c.append(calctime/steps)

for i in range(0, N):
  w_a = [];
  w_b = [];
  w_c = [];
  for j in range(0, CN):
    w_a.append(a[i*CN + j])
    w_b.append(b[i*CN + j])
    w_c.append(c[i*CN + j])
  wa_a = np.mean(w_a)
  wa_b = np.mean(w_b)
  wa_c = np.mean(w_c)
  print '{:8.0f}: {:8.0f} {:8.2f} {:8.2f} {:8.2f}                     {:8.2f}  {:8.2f}'.format(i + 1, steps, wa_c, wa_a, wa_c + wa_a, wa_b, wa_a + wa_b + wa_c)


print '--------------------------------------------------------------------------------------'
print ' Average: {:8.0f} {:8.2f} {:8.2f} {:8.2f}                     {:8.2f}  {:8.2f}'.format(steps, np.mean(c), np.mean(a), np.mean(c) + np.mean(a), np.mean(b), np.mean(a) + np.mean(b) + np.mean(c))
print '--------------------------------------------------------------------------------------'




