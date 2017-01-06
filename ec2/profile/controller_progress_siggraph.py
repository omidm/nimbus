#!/usr/bin/env python

import sys
import argparse
import os
import re
import numpy as np
import time


## Parse the command line arguments ##
parser = argparse.ArgumentParser(description='Process controller log file.')
parser.add_argument(
    "-i", "--input",
    dest="ifname",
    help="the input file to be processed, it should be the stdout of the controller")
parser.add_argument(
    "-d", "--dir_path",
    dest="dirpath",
    help="directory that holds the downloded logs, if provided overrides the -i option and works with the fixed structure of the folder")
parser.add_argument(
    "-s", "--start_index",
    dest="startindex",
    default=0,
    help="the first interation index to be included in the mean.")
parser.add_argument(
    "-e", "--end_index",
    dest="endindex",
    default=0,
    help="the last interation index to be included in the mean.")
parser.add_argument(
    "-p", "--plot",
    dest="plot",
    action="store_true",
    help="plot iteration length diagram")
parser.add_argument(
    "-v", "--verbose",
    dest="verbose",
    action="store_true",
    help="print per iteration stats as well")
parser.add_argument(
    "-a", "--assign_period",
    dest="assignperiod",
    action="store_true",
    help="iteration period is marked by assign, default is job done periods")

args = parser.parse_args()

SI = int(args.startindex)
EI = int(args.endindex)
assert(EI >= SI or EI == 0)

if (args.dirpath is not None):
  file_name = args.dirpath + '/logs/controller/stdout'
elif (args.ifname is not None):
  file_name = args.ifname
else:
  print "\nERROR: provide the input file or input folder!\n" 
  parser.print_help()
  exit(1)

PERIOD_MARK = "MARK_JOBDONE"
if (args.assignperiod):
  PERIOD_MARK = "MARK_ASSIGN"


f = open(file_name, 'r')

iteration_length     = []
version_manager      = []
add_complex_job      = []
data_manager_query   = []
data_manager_update  = []
template_instantiate = []
projection_time      = []
projection_state = 0
iter_num         = 0

for line in f:

  if projection_state == 1:
    if "projection_main" in line:
      projection_state = 2
  elif projection_state == 2:
    if "DYNAMIC" in line:
      projection_state = 3
      result = re.findall('DYNAMIC: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) .*', line)
      assert(len(result) > 0)
      start_projection_time = float(result[0])
  elif projection_state == 3:
    if "projection_loop_iteration_end" in line:
      projection_state = 4
  elif projection_state == 4:
    if "DYNAMIC" in line:
      projection_state = 0
      result = re.findall('DYNAMIC: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) .*', line)
      assert(len(result) > 0)
      projection_time[iter_num - 1] += (float(result[0]) - start_projection_time)

  if " main" in line:
    result = re.findall('.*: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) Picked .* main\.', line)
    if len(result) > 0:
      version_manager.append(0)
      add_complex_job.append(0)
      data_manager_query.append(0)
      data_manager_update.append(0)
      template_instantiate.append(0)
      last_time_stamp = float(result[0])
      continue

  if PERIOD_MARK in line:
    result = re.findall('.*: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) .*$', line)
    if len(result) > 0:
      iter_num += 1
      length = float(result[0]) - last_time_stamp
      iteration_length.append(length)
      version_manager.append(0)
      add_complex_job.append(0)
      data_manager_query.append(0)
      data_manager_update.append(0)
      template_instantiate.append(0)
      projection_time.append(0)
      projection_state = 1
      last_time_stamp = float(result[0])
      continue

  if "COMPLEX: VersionManager" in line:
    result = re.findall('COMPLEX: VersionManager: .* (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+)$', line)
    if len(result) > 0:
      version_manager[iter_num] += float(result[0])
      continue

  if "TEMPLATE: COMPLEX SPAWN" in line:
    result = re.findall('TEMPLATE: COMPLEX SPAWN: .* (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+)$', line)
    if len(result) > 0:
      add_complex_job[iter_num] += float(result[0])
      continue

  if "COMPLEX: QueryDataManager" in line:
    result = re.findall('COMPLEX: QueryDataManager: .* (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+)$', line)
    if len(result) > 0:
      data_manager_query[iter_num] += float(result[0])
      continue

  if "COMPLEX: Instantiate" in line:
    result = re.findall('COMPLEX: Instantiate: .* (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+)$', line)
    if len(result) > 0:
      template_instantiate[iter_num] += float(result[0])
      continue

  if "COMPLEX: UpdateDataManager" in line:
    result = re.findall('COMPLEX: UpdateDataManager: .* (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+)$', line)
    if len(result) > 0:
      data_manager_update[iter_num] += float(result[0])
      continue

version_manager.pop(iter_num)
add_complex_job.pop(iter_num)
data_manager_query.pop(iter_num)
data_manager_update.pop(iter_num)
template_instantiate.pop(iter_num)

assert(len(iteration_length) == iter_num)
assert(len(data_manager_update) == iter_num)

f.close()

if EI == 0:
  EI = iter_num

print "*** Iteration Number         : " + str(iter_num)
print "*** Considering iterations " + str(SI + 1) + " to " + str(EI) + " for the average stats:"
print "*** Average Iteration Length : " + str(np.mean(iteration_length[SI:EI]))

print '-------------------------------------------------------------------------------------------------------------'
print '       iter#       AC(VM)     DM Query    Temp Inst.   DM Update   (AC+DMQ+TI)   total    Proj.    Iter Leng'

if (args.verbose):
  print '-------------------------------------------------------------------------------------------------------------'
  for i in range(0, iter_num):
    print '    {:8.0f} {:8.3f}({:0.2f})   {:8.3f}     {:8.3f}     {:8.3f}     {:8.3f} {:8.3f} {:8.3f}   {:8.3f}'.format(
        i+1,
        add_complex_job[i],
        version_manager[i],
        data_manager_query[i],
        template_instantiate[i],
        data_manager_update[i],
        add_complex_job[i] + data_manager_query[i] +  template_instantiate[i],
        add_complex_job[i] + data_manager_query[i] +  template_instantiate[i] +data_manager_update[i],
        projection_time[i],
        iteration_length[i])


VM  = np.mean(version_manager[SI:EI])
AC  = np.mean(add_complex_job[SI:EI])
DMQ = np.mean(data_manager_query[SI:EI])
TIN = np.mean(template_instantiate[SI:EI])
DMU = np.mean(data_manager_update[SI:EI])
PT = np.mean(projection_time[SI:EI])
IL = np.mean(iteration_length[SI:EI])
print '-------------------------------------------------------------------------------------------------------------'
print 'Average: {:3.0f} {:8.3f}({:0.2f})   {:8.3f}     {:8.3f}     {:8.3f}     {:8.3f} {:8.3f} {:8.3f}   {:8.3f}'.format(
    EI - SI, AC, VM, DMQ, TIN, DMU, AC+DMQ+TIN, AC+DMQ+TIN+DMU, PT, IL) 
print '-------------------------------------------------------------------------------------------------------------'



if (not args.plot):
  exit(0)


import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# time_array = [t / 60 for t in time_array]
# 
# iter_num = range(0, len(time_array))
# 
# line = plt.plot(time_array, iter_num, '-*')
# plt.xlabel('Time (minute)')
# plt.ylabel('Iteration Number')
# # plt.axis([40, 160, 0, 0.03])
# plt.grid(True)
# plt.show()



line = plt.bar(range(1, iter_num+1), iteration_length)
plt.xlabel('Iteration Number')
plt.ylabel('Iteration Duration (second)')
# plt.axis([40, 160, 0, 0.03])
plt.grid(True)
plt.show()



