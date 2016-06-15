#!/usr/bin/python

# ./get_summary.py [directory_name] [number of workers]

import sys
import argparse
import numpy as np

def parse_line(line):
    items = line.split()
    return float(items[0]), float(items[3]), float(items[5]), float(items[7])

def sum(list):
  sum = 0
  for element in list:
    sum += element
  return sum

def print_bar(list, max_len):
  m = max(list)
  idx = 0;
  for element in list:
    idx += 1
    rod = "{:3d} : {:0.3f} *".format(idx, element)
    for i in range(0, int(max_len * element / m)):
      rod += "*"
    print rod

## Parse the command line arguments ##
parser = argparse.ArgumentParser(description='Process log files.')
parser.add_argument(
    "-d", "--dir_path",
    dest="dirpath",
    help="input dorectory path")
parser.add_argument(
    "-ti", "--truncate_index",
    dest="truncateindex",
    default=0,
    help="truncate index to ignore the initial iterations")
parser.add_argument(
    "-v", "--verbose",
    dest="verbose",
    action="store_true",
    help="print per iteration stats as well")
parser.add_argument(
    "-bl", "--bar_length",
    dest="barlength",
    default=50,
    help="the maximum length of the the printed bar per iteration that shows iteration length")

args = parser.parse_args()
d  = args.dirpath 
TI = int(args.truncateindex)
BL = int(args.barlength)





f = open('{}/controller_stats.txt'.format(d), 'r')

bytes_sent       = [];
bytes_received   = [];
overhead_time    = [];
iteration_length = []
last_time_stamp  = 0;
iter_num         = 0


for line in f:
    iter_num  += 1
    stamp, sent, received, overhead = parse_line(line)

    if (iter_num > 1):
      iteration_length.append(stamp - last_time_stamp)

    last_time_stamp = stamp
    overhead_time.append(overhead)
    bytes_sent.append(sent)
    bytes_received.append(received)

# iteration_length = [np.mean(iteration_length)] + iteration_length
iteration_length = [-1] + iteration_length
print iteration_length

if (args.verbose):
  print_bar(iteration_length, BL)


print "---------------------------------------------------------------"
print "Iteration number:  {:d}".format(iter_num)
print "** Truncated the first {:d} iteartions for the following stats:".format(TI)
print "---------------------------------------------------------------"

print "Average iteration duration | {:8.3f} (s)".format(np.mean(iteration_length[TI: iter_num]))
print "Average overhead           | {:8.3f} (s)".format(np.mean(overhead_time[TI: iter_num]))
print "Average bytes sent         | {:8.3f} (MB)".format(np.mean(bytes_sent[TI: iter_num]))
print "Average bytes received     | {:8.3f} (MB)".format(np.mean(bytes_received[TI: iter_num]))
print "---------------------------------------------------------------"

print "---------------------------------------------------------------"
print "Total bytes sent           | {:8.3f} (MB)".format(sum(bytes_sent))
print "Total bytes received       | {:8.3f} (MB)".format(sum(bytes_received))
print "---------------------------------------------------------------"

