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
        print drift
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
          print type(result[0][0])
      else:
        print "Wrong input in the file:" + file_name

  f.close()


for 




exit(0)

# print "Total " + args.tag + " time: " + str(time)

diff = []
for i in range (1, len(time) - 1):
  diff.append(float(time[i + 1] - time[i]))

diff.remove(max(diff))
diff.remove(min(diff))

print min(diff)
print float(sum(diff)) / float(len(diff))
print max(diff)


import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# the histogram of the data
n, bins, patches = plt.hist(diff, 50, normed=1, facecolor='green', alpha=0.75)

plt.xlabel('Projection Period')
plt.ylabel('Probability')
# plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.show()



