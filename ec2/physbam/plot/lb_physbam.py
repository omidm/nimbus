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
base_time = 0
first = True

# How to count number od iterations from std::out of controller? 
# either count loop_iteration_part_two jobs: 
#     cat log | grep complex | grep "projection_loop_iteration_end\." -c
#                      +
#     cat omid | grep Picked | grep " loop_iteration_part_two\."
#
# or count loop_iteration jobs (loop_frame is not templetized):
#    cat omid | grep complex | grep "loop_iteration_part_two\." -c
#                     +
#    cat omid | grep Picked | grep " loop_iteration\." -c
#

for line in content:
  # result = re.findall('(\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) : ' + args.tag, line)
  result = re.findall('DT: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+)$', line)
  if len(result) > 0:
    time.append(decimal.Decimal(result[0]))
    if first:
      base_time = decimal.Decimal(result[0])
      first = False

f.close()

print "Time: "
print time
print "Base: "
print base_time

for i in range (0, len(time)):
  time[i] = time[i] - base_time

diff = []
for i in range (0, len(time) - 1):
  diff.append(time[i + 1] - time[i])

print "Time: "
print time
print "Average duration: " + str(numpy.mean(diff))
print "Iteration Number: " + str(len(diff))

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

time = [t / 60 for t in time]

iter_num = range(0, len(time))

line = plt.plot(time, iter_num, '-*')

plt.xlabel('Time (minute)')
plt.ylabel('Iteration Number')
# plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.show()


iter_num = range(0, len(time) - 1)

line = plt.bar(iter_num, diff)

plt.xlabel('Iteration Number')
plt.ylabel('Iteration Duration (second)')
# plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.show()



