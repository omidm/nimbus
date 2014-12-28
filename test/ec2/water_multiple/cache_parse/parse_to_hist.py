#!/usr/bin/env python 

import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser
import sys

###############################################################################
#                               PARSER OPTIONS                                #
###############################################################################

parser = OptionParser()
parser.add_option('-i', '--in', dest='input', help='file containing data to parse')
parser.add_option('-o', '--out', dest='output', help='file to store results')
parser.add_option('-p', '--plot', dest='plt_type', help='plot pdf or cdf')

(options, args) = parser.parse_args()

###############################################################################
#                                    PARSE                                    #
###############################################################################

with open(options.input) as data:
  num_lines = len(data.readlines())

data = open(options.input)

# store start times to calculate delta, and delta times to plot histogram
start_times = dict()
delta_times = dict()

# parse each line, make state transitions and store times
num_line = 1
percent_check = 2
percent_inc = 2
sys.stdout.write("Parsed file % :   0")
sys.stdout.flush()
for line in data:
  words = line.split(";")
  thread = words[0].strip()
  status = words[1].strip()
  event = words[2].strip()
  time = float(words[3].strip())
  if status == 'start':
    if thread not in start_times.keys():
      start_times[thread] = dict()
    start_times[thread][event] = time
  else:
    if status == 'end':
      if event not in delta_times.keys():
        delta_times[event] = []
      delta_times[event] += [time - start_times[thread][event]]
  num_line += 1
  if num_line * 100 / num_lines == percent_check:
    sys.stdout.write("\b\b\b%3i" % percent_check)
    sys.stdout.flush()
    percent_check += percent_inc
sys.stdout.write("\n")
sys.stdout.flush()

# plot histograms
print("Plotting %s" % options.plt_type)
num_events = len(delta_times.keys())
events = delta_times.keys()
events.sort()
subplot = 1
for event in events:
  plt.subplot((num_events+1)/2, 2, subplot)
  plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  if options.plt_type == 'pdf':
    plt.hist(delta_times[event], np.linspace(0, 3e-5, 500), \
            color='black', alpha=0.15)
    plt.ylabel('Frequency')
  elif options.plt_type == 'cdf':
    plt.hist(delta_times[event], np.linspace(0, 3e-5, 500), \
            normed=True, cumulative=True, histtype='step', \
            color='black', alpha=0.1)
    plt.ylabel('CDF')
  plt.title(event)
  plt.xlabel('Time (s)')
  subplot += 1
plt.subplots_adjust(hspace=1.25)
plt.savefig(options.output)
