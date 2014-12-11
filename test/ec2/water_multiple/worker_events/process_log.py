#!/usr/bin/env python

import decimal
import argparse
import os
import os.path
import re
import numpy
import time

expr = re.compile("[^0-9\.\-]+")
def get_key(line):
    ls = re.split(expr, line)
    if ls[0]:
        return float(ls[0])
    else:
        return float(ls[1])


## Parse the command line arguments ##
parser = argparse.ArgumentParser(description='Process log files.')
# parser.add_argument(
#     "-i", "--input",
#     dest="input_base_name",
#     default="log",
#     help="base input file name that follows with -<worker number>")
parser.add_argument(
    "-d", "--directory",
    dest="dir",
    default=".",
    help="directory to find the input files")
parser.add_argument(
    "-od", "--outdirectory",
    dest="odir",
    default=".",
    help="directory to write the outout files")
parser.add_argument(
    "-wn", "--workernum",
    dest="worker_num",
    default=1,
    type=int,
    help="number of workers to process")
parser.add_argument(
    "-cn", "--corenum",
    dest="core_num",
    default=1,
    type=int,
    help="number of cores per worker")
parser.add_argument(
    "-o", "--output",
    dest="output_base_name",
    default="data",
    help="base output file name that follows with -<worker number>")

args = parser.parse_args()

CN = args.core_num

for i in range(1, args.worker_num + 1):
  print 'Processing logs for worker ' + str(i) + '...'

  fe_log_file = args.dir + '/' + str(i) + '_event_fe.txt'
  be_log_file = args.dir + '/' + str(i) + '_event_be.txt'
  data_file =  args.odir + '/' + args.output_base_name + '-' + str(i)

  print 'Opening file ' + fe_log_file
  fe_log = open(fe_log_file, 'r')
  print 'Opening file ' + be_log_file
  be_log = open(be_log_file, 'r')

  start_time = time.time()
  print 'Reading in files ...'
  fe_lines = fe_log.readlines()
  for i in range(1, 1000):
    fe_lines.pop()
  be_lines = be_log.readlines()
  for i in range(1, 1000):
    be_lines.pop()

  lines = numpy.concatenate((fe_lines, be_lines), axis=0)

  end_time = time.time()
  print start_time - end_time

  start_time = time.time()
  print 'Sorting lines ...'
  sorted_lines = sorted(lines, key=get_key)
  end_time = time.time()
  print start_time - end_time


  total_lines = len(sorted_lines)
  print "Parsing %d lines ..." % total_lines

  ### Initialize data to get from parser ###
  copy_time    = 0
  compute_time = 0
  running_time = 0
  idle_time    = 0
  blocked_time = 0

  blocked_jobs = 0
  ready_jobs   = 0
  running_jobs = 0

  init_time_stamp = float(re.findall('(\d+\.\d+) .*', sorted_lines[0])[0])
  last_time_stamp = init_time_stamp

  start_time = time.time()
  num = 0;
  for line in sorted_lines:
      num += 1
      if num == total_lines/4:
          print "Parsed 25% ..."
          end_time = time.time()
          print start_time - end_time
          start_time = end_time
      if num == total_lines/2:
          print "Parsed 50% ..."
          end_time = time.time()
          print start_time - end_time
          start_time = end_time
      if num == 3*total_lines/4:
          print "Parsed 75% ..."
          end_time = time.time()
          print start_time - end_time
          start_time = end_time


      if "JOB_DONE" in line:
        continue

      parsed = []
      if " dispatch_job(job_done) " in line:
        parsed = re.findall('(\d+\.\d+) .* (\d+) \d+', line)
      else:
        parsed = re.findall('(\d+\.\d+) .* (\d+)', line)

      if len(parsed) > 0:
        time_stamp = float(parsed[0][0])
        job_num    = int(parsed[0][1])
      else:
        print "ERROR: bad input line: " + line
        continue



      elapsed = time_stamp - last_time_stamp
      running_time += (elapsed * running_jobs)
      if "Copy" in line:
        copy_time += (elapsed * running_jobs)
      else:
        compute_time += (elapsed * running_jobs)
        
      if running_jobs < CN:
        blocked_time += (elapsed * min(blocked_jobs, CN - running_jobs))

      if (running_jobs + blocked_jobs) < CN:
        idle_time += (elapsed * (CN - running_jobs - blocked_jobs))


      if " recv_job " in line:
        blocked_jobs += 1

      elif " dispatch_job(" in line:
        blocked_jobs -= 1
        ready_jobs += 1

      elif " r " in line:
        ready_jobs -= 1
        running_jobs += 1

      elif " f " in line:
        running_jobs -= 1


      last_time_stamp = time_stamp


  print "Parsed 100% ..."
  print blocked_jobs
  print ready_jobs
  print running_jobs

  
  print "Opening %s ..." % data_file
  data = open(data_file, 'w')
  data.write("Total Wall Time %2.6f\n" % (last_time_stamp - init_time_stamp))
  data.write("Total Cycle Time %2.6f\n" % (blocked_time + running_time + idle_time))
  data.write("Blocked Time %2.6f\n" % blocked_time)
  data.write("Running Time %2.6f\n" % running_time)
  data.write("    Copy Time %2.6f\n" % copy_time)
  data.write("    Compute Time %2.6f\n" % compute_time)
  data.write("Idle Time %2.6f\n" % idle_time)

  data.close()
  fe_log.close()
  be_log.close()


