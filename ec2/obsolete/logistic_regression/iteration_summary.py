#!/usr/bin/python

# ./get_summary.py [directory_name] [number of workers]

import sys
import argparse

def parse_line(line):
    items = line.split()
    pexec = float('nan')
    dxl   = float('nan')
    if len(items) == 8:
      pexec = float(items[7])
    elif len(items) == 10:
      pexec = float(items[7])
      dxl   = float(items[9])
    elif len(items) == 12:
      pexec = float(items[11])
      dxl   = float(items[9])
    return float(items[1]),  float(items[3]), float(items[5]), pexec, dxl


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
    "-ti", "--truncate_index",
    dest="truncateindex",
    default=0,
    help="truncate index to ignore the initial iterations")
parser.add_argument(
    "-v", "--verbose",
    dest="collapse",
    action="store_true",
    help="print per iteration stats as well")


args = parser.parse_args()

d  = args.dirpath 
N  = int(args.workernum) 
CN = int(args.corenum)
TI = int(args.truncateindex)





total_dxl_sum    = 0;
total_parent_sum = 0;
total_run_sum    = 0;
total_block_sum  = 0;
total_idle_sum   = 0;
total_total_sum  = 0;
iter_nums = []
    
print '--------------------------------------------------------------------------------------'
print '   worker    steps   parent  dx_lock      run     block     idle    synch     total'
print '--------------------------------------------------------------------------------------'

for rank in range(1, N+1):
  f = open('{}/{}_main_timers.txt'.format(d, rank), 'r')

  iter_num   = 0
  dxl_sum    = 0;
  parent_sum = 0;
  run_sum    = 0;
  block_sum  = 0;
  idle_sum   = 0;
  total_sum  = 0;

  idx = 0
  for line in f:
      idx += 1
      run, block, idle, parent, dxl = parse_line(line)
      if (args.collapse):
        print '          {:8.0f} {:8.3f} {:8.3f} {:8.3f}  {:8.3f} {:8.3f} {:8.3f}  {:8.3f}'.format(
            idx, parent, dxl, run/CN, block/CN, idle/CN, (block+idle)/CN, (idle+block+run)/CN)
      if idx <= TI:
        continue
      iter_num  += 1
      dxl_sum    += dxl
      parent_sum += parent
      run_sum    += run
      block_sum  += block
      idle_sum   += idle
      total_sum  += run + block + idle

  iter_nums.append(iter_num)
  total_dxl_sum    += dxl_sum/iter_num
  total_parent_sum += parent_sum/iter_num
  total_run_sum    += run_sum/iter_num
  total_block_sum  += block_sum/iter_num
  total_idle_sum   += idle_sum/iter_num
  total_total_sum  += total_sum/iter_num

  print '{:8.0f}: {:8.0f} {:8.3f} {:8.3f} {:8.3f}  {:8.3f} {:8.3f} {:8.3f}  {:8.3f}'.format(
            rank, iter_num, parent_sum/iter_num, dxl_sum/iter_num, run_sum/CN/iter_num, block_sum/CN/iter_num, idle_sum/CN/iter_num, (block_sum+idle_sum)/CN/iter_num, total_sum/CN/iter_num)

print '--------------------------------------------------------------------------------------'
print ' Average: {:8.0f} {:8.3f} {:8.3f} {:8.3f}  {:8.3f} {:8.3f} {:8.3f}  {:8.3f}'.format(
    iter_nums[0], total_parent_sum, total_dxl_sum/N, total_run_sum/N/CN, total_block_sum/N/CN, total_idle_sum/N/CN, (total_block_sum+total_idle_sum)/N/CN, total_total_sum/N/CN)
print '--------------------------------------------------------------------------------------'
