#!/usr/bin/python

import numpy as np
import sys
import math
import argparse
import pylab
import matplotlib.pyplot as plt
from matplotlib import rc

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

task_length = []

for rank in range(1, CN*N + 1):
    packtime   = 0
    synctime   = 0
    calcsubstep    = 0
    init_phase = True
    f = open('{}mpi{}.log'.format(d, rank))

    steps = 0
    for line in f:
        if   'Pack' in line \
          or 'Prepare' in line \
          or 'Unpack' in line \
          or 'Wait' in line \
          or 'End' in line \
          or 'ImplicitPart' in line \
          or 'ProjectionGlobal' in line \
          or 'CalculateSubstep' in line \
          or 'CalculateFrame' in line \
          or 'Simulation ' in line:
            continue;
        task_length.append(float(line.split()[1]));




print '---------------------------PhysBAM Task Length Distribution---------------------------'
print ' {:20s} : {:8.0f} '.format('sample number', len(task_length))
print ' {:20s} : {:8.7f} '.format('maximum', np.max(task_length))
print ' {:20s} : {:8.7f} '.format('minmum', np.min(task_length))
print ' {:20s} : {:8.7f} '.format('average', np.mean(task_length))
print ' {:20s} : {:8.7f} '.format('5th  percentile', np.percentile(task_length, 5))
print ' {:20s} : {:8.7f} '.format('10th percentile', np.percentile(task_length, 10))
print ' {:20s} : {:8.7f} '.format('20th percentile', np.percentile(task_length, 20))
print ' {:20s} : {:8.7f} '.format('30th percentile', np.percentile(task_length, 30))
print ' {:20s} : {:8.7f} '.format('40th percentile', np.percentile(task_length, 40))
print ' {:20s} : {:8.7f} '.format('50th percentile', np.percentile(task_length, 50))
print ' {:20s} : {:8.7f} '.format('60th percentile', np.percentile(task_length, 60))
print ' {:20s} : {:8.7f} '.format('70th percentile', np.percentile(task_length, 70))
print ' {:20s} : {:8.7f} '.format('80th percentile', np.percentile(task_length, 80))
print ' {:20s} : {:8.7f} '.format('90th percentile', np.percentile(task_length, 90))
print ' {:20s} : {:8.7f} '.format('95th percentile', np.percentile(task_length, 95))
print ' {:20s} : {:8.7f} '.format('98th percentile', np.percentile(task_length, 98))
print ' {:20s} : {:8.7f} '.format('99th percentile', np.percentile(task_length, 99))
print '--------------------------------------------------------------------------------------'


rc('font', size=24)

num_bins = 400
cut_percentile = 100

# bins = np.linspace(np.min(task_length), np.percentile(task_length, cut_percentile), num=num_bins)
bins = np.logspace(math.log(np.min(task_length), 10), math.log(np.percentile(task_length, cut_percentile), 10), num=num_bins)
weights = np.ones_like(task_length)/float(len(task_length))


pylab.figure()
pylab.hist(task_length, linewidth=3, bins=bins, weights=weights, normed=False, histtype='step', stacked=True, fill=False, cumulative=True)

pylab.xlabel('Task length (seconds)', fontsize='x-large')
pylab.xticks(fontsize='x-large')
# pylab.xticks(np.linspace(0,.070,8), fontsize='large')
# pylab.xlim([0, .071])
pylab.ylabel('CDF', fontsize='x-large')
pylab.yticks(fontsize='large')
# pylab.yticks(np.linspace(0,0.4,5), fontsize='large')
pylab.xscale('log')
pylab.xlim([0, 15])

pylab.show()



# counts, bin_edges = np.histogram(task_length, bins=num_bins, normed=True)
# cdf = np.cumsum(counts)
# 
# plt.xscale('log')
# 
# plt.plot(bin_edges[1:], cdf)
# plt.show()














