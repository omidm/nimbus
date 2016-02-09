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
    "-d", "--directory",
    dest="idir",
    default=".",
    help="directory to find the input file")
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
parser.add_argument(
    "-t", "--tag",
    dest="tag",
    default="assign",
    help="tag to sum up the time for")


args = parser.parse_args()

file_name = args.idir + '/' + args.ifname
print 'Opening the file ' + file_name

f = open(file_name, 'r')
content = f.readlines()

# time = 0
time = [];

for line in content:
  # result = re.findall('(\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) : ' + args.tag, line)
  result = re.findall('INFO: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) id:.*', line)
  if len(result) > 0:
    # time += decimal.Decimal(result[0])
    time.append(decimal.Decimal(result[0]))

f.close()

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



