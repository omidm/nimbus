#!/usr/bin/env python

import sys
import argparse
import os
import re
import numpy as np
import decimal
import time


## Parse the command line arguments ##
parser = argparse.ArgumentParser(description='Process log files.')
parser.add_argument(
    "-d", "--dir_path",
    dest="dirpath",
    default=".",
    help="input dorectory path")
parser.add_argument(
    "-i", "--input",
    dest="ifname",
    default="ec2_log.txt",
    help="input file name")
parser.add_argument(
    "-ti", "--truncate_index",
    dest="truncateindex",
    default=0,
    help="truncate index to start the mean from there")
parser.add_argument(
    "-p", "--plot",
    dest="plot",
    action="store_true",
    help="plot iteration length plots")
parser.add_argument(
    "-v", "--verbose",
    dest="collapse",
    action="store_true",
    help="print per iteration stats as well")

args = parser.parse_args()

TI = int(args.truncateindex)

file_name = args.dirpath + '/' + args.ifname
# print 'Opening the file ' + file_name

f = open(file_name, 'r')

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

iteration_length     = []
version_manager      = []
add_complex_job      = []
data_manager_query   = []
data_manager_update  = []
template_instantiate = []
iter_num         = 0

for line in f:
  # result = re.findall('(\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) : ' + args.tag, line)
  if " main" in line:
    result = re.findall('.*: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) Picked .* main\.', line)
    if len(result) > 0:
      version_manager.append(0)
      add_complex_job.append(0)
      data_manager_query.append(0)
      data_manager_update.append(0)
      template_instantiate.append(0)
      last_time_stamp = decimal.Decimal(result[0])
      continue
  if "__MARK_STAT" in line:
    result = re.findall('.*: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) complex .* __MARK_STAT.*\.', line)
    if len(result) > 0:
      iter_num += 1
      length = decimal.Decimal(result[0]) - last_time_stamp
      iteration_length.append(length)
      version_manager.append(0)
      add_complex_job.append(0)
      data_manager_query.append(0)
      data_manager_update.append(0)
      template_instantiate.append(0)
      last_time_stamp = decimal.Decimal(result[0])
      continue
    result = re.findall('.*: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) Picked .* __MARK_STAT.*\.', line)
    if len(result) > 0:
      iter_num += 1
      length = decimal.Decimal(result[0]) - last_time_stamp
      iteration_length.append(length)
      version_manager.append(0)
      add_complex_job.append(0)
      data_manager_query.append(0)
      data_manager_update.append(0)
      template_instantiate.append(0)
      last_time_stamp = decimal.Decimal(result[0])
      continue

  #############################################################################################################
  # This part is obsolete, just kept it so that old logs could be processed
  if " loop_iteration" in line:
    result = re.findall('.*: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) complex .* loop_iteration_part_two\.', line)
    if len(result) > 0:
      iter_num += 1
      length = decimal.Decimal(result[0]) - last_time_stamp
      iteration_length.append(length)
      version_manager.append(0)
      add_complex_job.append(0)
      data_manager_query.append(0)
      data_manager_update.append(0)
      template_instantiate.append(0)
      last_time_stamp = decimal.Decimal(result[0])
      continue
    result = re.findall('.*: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) Picked .* loop_iteration\.', line)
    if len(result) > 0:
      iter_num += 1
      length = decimal.Decimal(result[0]) - last_time_stamp
      iteration_length.append(length)
      version_manager.append(0)
      add_complex_job.append(0)
      data_manager_query.append(0)
      data_manager_update.append(0)
      template_instantiate.append(0)
      last_time_stamp = decimal.Decimal(result[0])
      continue
  #############################################################################################################

  if "COMPLEX: VersionManager" in line:
    result = re.findall('COMPLEX: VersionManager: .* (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+)$', line)
    if len(result) > 0:
      version_manager[iter_num] += decimal.Decimal(result[0])
      continue
  if "TEMPLATE: COMPLEX SPAWN" in line:
    result = re.findall('TEMPLATE: COMPLEX SPAWN: .* (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+)$', line)
    if len(result) > 0:
      add_complex_job[iter_num] += decimal.Decimal(result[0])
      continue
  if "COMPLEX: QueryDataManager" in line:
    result = re.findall('COMPLEX: QueryDataManager: .* (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+)$', line)
    if len(result) > 0:
      data_manager_query[iter_num] += decimal.Decimal(result[0])
      continue
  if "COMPLEX: Instantiate" in line:
    result = re.findall('COMPLEX: Instantiate: .* (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+)$', line)
    if len(result) > 0:
      template_instantiate[iter_num] += decimal.Decimal(result[0])
      continue
  if "COMPLEX: UpdateDataManager" in line:
    result = re.findall('COMPLEX: UpdateDataManager: .* (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+)$', line)
    if len(result) > 0:
      data_manager_update[iter_num] += decimal.Decimal(result[0])
      continue

version_manager.pop(iter_num)
add_complex_job.pop(iter_num)
data_manager_query.pop(iter_num)
data_manager_update.pop(iter_num)
template_instantiate.pop(iter_num)

assert(len(iteration_length) == iter_num)
assert(len(data_manager_update) == iter_num)

f.close()


print "*** Iteration Number         : " + str(iter_num)
print "*** Truncated the first " + str(TI) + " iteartions for the average stats:"
print "*** Average Iteration Length : " + str(np.mean(iteration_length[TI:iter_num-1]))

print '----------------------------------------------------------------------------------------------------'
print '       iter#       AC(VM)     DM Query    Temp Inst.   DM Update   (AC+DMQ+TI)   total   Iter Leng'

if (args.collapse):
  print '----------------------------------------------------------------------------------------------------'
  for i in range(0, iter_num):
    print '    {:8.0f} {:8.2f}({:0.2f})   {:8.2f}     {:8.2f}     {:8.2f}     {:8.2f} {:8.2f}    {:8.2f}'.format(
        i+1,
        add_complex_job[i],
        version_manager[i],
        data_manager_query[i],
        template_instantiate[i],
        data_manager_update[i],
        add_complex_job[i] + data_manager_query[i] +  template_instantiate[i],
        add_complex_job[i] + data_manager_query[i] +  template_instantiate[i] +data_manager_update[i],
        iteration_length[i])


VM  = np.mean(version_manager[TI:iter_num])
AC  = np.mean(add_complex_job[TI:iter_num])
DMQ = np.mean(data_manager_query[TI:iter_num])
TIN = np.mean(template_instantiate[TI:iter_num])
DMU = np.mean(data_manager_update[TI:iter_num])
IL = np.mean(iteration_length[TI:iter_num])
print '----------------------------------------------------------------------------------------------------'
print 'Average: {:3.0f} {:8.2f}({:0.2f})   {:8.2f}     {:8.2f}     {:8.2f}     {:8.2f} {:8.2f}    {:8.2f}'.format(
    iter_num - TI, AC, VM, DMQ, TIN, DMU, AC+DMQ+TIN, AC+DMQ+TIN+DMU, IL) 
print '----------------------------------------------------------------------------------------------------'




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



