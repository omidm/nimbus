#!/usr/bin/python

# ./get_summary.py [directory_name] [number of workers]

import sys

def parse_line(line):
    items = line.split()
    return float(items[1]),  float(items[3]), float(items[5])

total_run_sum   = 0;
total_block_sum = 0;
total_idle_sum  = 0;
total_total_sum = 0;
iter_nums = []
    
print '   worker    steps      run    block     idle    total'

N = int(sys.argv[2])
for rank in range(1, N+1):
  f = open('{}/{}_main_timers.txt'.format(sys.argv[1], rank), 'r')

  iter_num  = 0
  run_sum   = 0;
  block_sum = 0;
  idle_sum  = 0;
  total_sum = 0;

  idx = 0
  for line in f:
      idx += 1
      if idx == 1:
        continue
      run, block, idle = parse_line(line)
      iter_num  += 1
      run_sum   += run
      block_sum += block
      idle_sum  += idle
      total_sum  += run + block + idle

  iter_nums.append(iter_num)
  total_run_sum   += run_sum
  total_block_sum += block_sum
  total_idle_sum  += idle_sum
  total_total_sum  += total_sum

  print '{:8.0f}: {:8.0f} {:8.2f} {:8.2f} {:8.2f} {:8.2f}'.format(
            rank, iter_num, run_sum/8, block_sum/8, idle_sum/8, total_sum/8)
print '---------------------------------------------------------'
print ' Average: {:8.0f} {:8.2f} {:8.2f} {:8.2f} {:8.2f}'.format(iter_nums[0], total_run_sum/N/8, total_block_sum/N/8, total_idle_sum/N/8, total_total_sum/N/8)
