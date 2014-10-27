#!/usr/bin/env python

import sys
import argparse
import os
import re
import numpy
import decimal


## Parse the command line arguments ##
parser = argparse.ArgumentParser(description='Process log files.')
parser.add_argument(
    "-i", "--input",
    dest="ifname",
    default="log",
    help="input file name")
parser.add_argument(
    "-od", "--outdirectory",
    dest="odir",
    default=".",
    help="directory to dump the output file")
parser.add_argument(
    "-o", "--output",
    dest="ofname",
    default="data",
    help="output file name")

args = parser.parse_args()

FN = 2
throughput = 1.008 # Gbps
rtt = 0.7e-3


STime = {}
RTime = {}
SSize = {}
RSize = {}

for n in range (1, FN + 1):

  file_name = args.ifname + str(n)
  print 'Opening the file ' + file_name
  
  f = open(file_name, 'r')
  content = f.readlines()

  Drift = []
  drift = 0

  for line in content:
    if 'D' in line:
      result = re.findall('D (-*\d+|-*\d+\.\d+|-*\d+e-\d+|-*\d+\.\d+e-\d+)$', line)
      if len(result) == 1:
        Drift.append(float(decimal.Decimal(result[0])))
        drift = numpy.mean(Drift)
      else:
        print "Wrong input in the file:" + file_name

    else:
      result = re.findall('(\w) (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) j:\s+(\d+) s:\s+(\d+)$', line)
      if len(result) == 1:
        if result[0][0] == 'S':
          STime[result[0][2]] = float(decimal.Decimal(result[0][1])) - drift
          SSize[result[0][2]] = float(decimal.Decimal(result[0][3]))
        elif result[0][0] == 'R':
          RTime[result[0][2]] = float(decimal.Decimal(result[0][1])) - drift
          RSize[result[0][2]] = float(decimal.Decimal(result[0][3]))
        else:
          print "Wrong input in the file:" + file_name
      else:
        print "Wrong input in the file:" + file_name

  f.close()


Speed = {}
Elapsed = {}
total_count = 0
negative_count = 0

for j in RTime:
  total_count += 1
  elapsed = RTime[j] - STime[j];
  if elapsed < 0:
    negative_count += 1
    elapsed = rtt
  assert(RSize[j] == SSize[j]);
  size = RSize[j]
  if Elapsed.has_key(size):
    Elapsed[size].append(elapsed)
    Speed[size].append(size * 8 / elapsed)
  else:
    Elapsed[size] = [elapsed]
    Speed[size] = [size * 8 / elapsed]

print "Different Size count: " + str(len(Elapsed))
print "Negative elapsed time count: " + str(negative_count)
print "Total elapsed time count: " + str(total_count)


SpeedMean = {}
SpeedError = {}
ElapsedMean = {}
ElapsedError = {}
Size = []

for s in Elapsed:
  speed_mean = numpy.mean(Speed[s])
  speed_error = numpy.std(Speed[s])
  elapsed_mean = numpy.mean(Elapsed[s])
  elapsed_error = numpy.std(Elapsed[s])
  if (speed_error > speed_mean):
    speed_error = speed_mean
  if (elapsed_error > elapsed_mean):
    elapsed_error = elapsed_mean
  Size.append(s)
  SpeedMean[s] = speed_mean
  SpeedError[s] = speed_error
  ElapsedMean[s] = elapsed_mean
  ElapsedError[s] = elapsed_error
  # print "size: " + str(s) + " spead mean: " + str(mean) + " speed std: " + str(error)

S = []
SC = []
SM = []
SE = []
EM = []
EE = []

Size.sort()
MainSize = Size[:]
  
for s in Size:
  S.append(s / 1000)
  # SM.append(SpeedMean[s])
  # SE.append(SpeedError[s])
  SM.append(s * 8 / ElapsedMean[s])
  SE.append(0)
  EM.append(ElapsedMean[s])
  EE.append(ElapsedError[s])
  SC.append(len(Elapsed[s]))

# bad sample!  
# EM[0] = rtt
# EE[0] = rtt / 2

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

font = FontProperties(family='sans-serif', weight='bold', size=12)


fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True)


newSC = map(lambda x: x/float(sum(SC)), SC)
ax2.bar(S, newSC, 1, color='g', hatch='', bottom=0, edgecolor="g")
ax2.set_ylabel('Distribution', family='sans-serif', size=12, weight='bold')


ax1.axhline(throughput, color='r', ls='--', linewidth=4)
tag = 'iperf throughput: %2.2f Gbps' % throughput
ax1.annotate(tag, xy=(S[len(S) / 2 + 2], throughput), fontproperties=font)

newSM = map(lambda x: x/1e9, SM)
newSE = map(lambda x: x/1e9, SE)
ax1.errorbar(S, newSM, yerr=newSE, fmt='-o')
ax1.set_ylabel('Transmission Speed (Gbps)', family='sans-serif', size=12, weight='bold')
ax1.grid(True)

ax0.axhline(rtt, color='r', ls='--', linewidth=4)
tag = 'ping RTT: %2.2f ms' % (rtt / 1e-3)
ax0.annotate(tag, xy=(S[len(S) / 2 + 2], rtt), fontproperties=font)

newEM = map(lambda x: x/1e-3, EM)
newEE = map(lambda x: x/1e-3, EE)
ax0.errorbar(S, newEM, yerr=newEE, fmt='-o')
ax0.set_ylabel('Transmission Time (ms)', family='sans-serif', size=12, weight='bold')
ax0.grid(True)

ax0.set_title('Analysis of Data Exchange for Water Simulation \n 512 cube, 64 application partitions, c3.2xlarge instances', family='sans-serif', size=12, weight='bold')
plt.xlabel('Data Size (KB)', family='sans-serif', size=12, weight='bold')
plt.xscale('linear')
plt.yscale('linear')
plt.rc('font', family='sans-serif', size=12, weight='bold')
# plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.show()



