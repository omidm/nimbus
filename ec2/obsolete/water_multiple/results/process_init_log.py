#!/usr/bin/env python

import decimal
import argparse
import os
import os.path
import re

## Parse the command line arguments ##
parser = argparse.ArgumentParser(description='Process log files.')
parser.add_argument(
    "-i", "--input",
    dest="input_base_name",
    default="log",
    help="base input file name that follows with -<worker number>")
parser.add_argument(
    "-d", "--directory",
    dest="dir",
    default=".",
    help="directory to find the input files")
parser.add_argument(
    "-od", "--outdirectory",
    dest="odir",
    default=".",
    help="directory to find the input files")
parser.add_argument(
    "-wn", "--workernum",
    dest="worker_num",
    default=1,
    type=int,
    help="number of workers to process")
parser.add_argument(
    "-o", "--output",
    dest="output_base_name",
    default="data",
    help="base output file name that follows with -<worker number>")

args = parser.parse_args()

for i in range(1, args.worker_num + 1):
  print 'Processing the log for worker ' + str(i) + '...'

  log_file = args.dir + '/' + args.input_base_name + '-' + str(i)
  data_file =  args.odir + '/' + args.output_base_name + '-' + str(i)

  print 'Opening file ' + log_file
  log_temp = open(log_file, 'r')
  log = open(log_file, 'r')
  
  total_lines = sum(1 for line in log_temp)
  print "Parsing %d lines ..." % total_lines

  ### Initialize data to get from parser ###
  init_total  = 0
  init_start  = 0
  init_end    = 0

  for num, line in enumerate(log):
      x =  re.findall('(\d+\.\d+$|\d+e-\d+$|\d+$)', line)
      if len(x) > 0:
          xnum = decimal.Decimal(x[0])
          if "start" in line:
              init_start = xnum
          if "end" in line:
              init_end = xnum
              init_total += (init_end - init_start) 
      if num == total_lines/4:
          print "Parsed 25% ..."
      if num == total_lines/2:
          print "Parsed 50% ..."
      if num == 3*total_lines/4:
          print "Parsed 75% ..."

  print "init time for worker %d: %2.2f " % (i, init_total)
  print "Opening %s ..." % data_file
  data = open(data_file, 'a')
  data.write("Initialize %2.2f\n" % init_total)

  data.close()
  log.close()
  log_temp.close()


