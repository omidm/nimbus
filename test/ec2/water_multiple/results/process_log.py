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
  worker_total  = 0
  compute_non_sterile_total = 0
  compute_sterile_total = 0
  copy_total    = 0
  projection_total = 0

  tnew_s_total_compute_non_sterile  = 0
  tnew_f_total_compute_non_sterile  = 0
  tnew_p_total_compute_non_sterile  = 0
  told_s_total_compute_non_sterile  = 0
  told_f_total_compute_non_sterile  = 0
  told_p_total_compute_non_sterile  = 0
  tnew_s_total_compute_sterile  = 0
  tnew_f_total_compute_sterile  = 0
  tnew_p_total_compute_sterile  = 0
  told_s_total_compute_sterile  = 0
  told_f_total_compute_sterile  = 0
  told_p_total_compute_sterile  = 0
  tnew_s_total_copy  = 0
  tnew_f_total_copy  = 0
  tnew_p_total_copy  = 0
  told_s_total_copy  = 0
  told_f_total_copy  = 0
  told_p_total_copy  = 0

  worker_start  = 0
  compute_non_sterile_start = 0
  compute_sterile_start = 0
  copy_start    = 0
  projection_start = 0

  tnew_s_start  = 0
  tnew_f_start  = 0
  tnew_p_start  = 0
  told_s_start  = 0
  told_f_start  = 0
  told_p_start  = 0

  worker_end    = 0
  compute_non_sterile_end   = 0
  compute_sterile_end   = 0
  copy_end      = 0
  projection_end = 0

  tnew_s_end    = 0
  tnew_f_end    = 0
  tnew_p_end    = 0
  told_s_end    = 0
  told_f_end    = 0
  told_p_end    = 0

  unknown_trans = 0
  latest_seen   = 0
  projection_trans = 0
  
  compute_non_sterile = False
  compute_sterile = False
  copy = False
  projection = False

  for num, line in enumerate(log):
      x =  re.findall('(\d+\.\d+$|\d+e-\d+$|\d+$)', line)
      if len(x) > 0:
          xnum = decimal.Decimal(x[0])
          latest_seen = xnum 
          if "Worker starts" in line:
              worker_start = xnum
          if "App compute job start" in line:
              if 'loop' in line:
                  compute_non_sterile_start = xnum
                  compute_non_sterile = True
              else:
                  compute_sterile_start = xnum
                  compute_sterile = True
                  if 'projection' in line:
                      projection_start = xnum
                      projection = True
          if "App copy job start" in line:
              copy_start = xnum
              copy = True
          if "Read Scalar Array (New Translator) start" in line or \
              "Write Scalar Array (New Translator) start" in line:
              tnew_s_start = xnum
          if "Read Face Array (New Translator) start" in line or \
              "Write Face Array (New Translator) start" in line:
              tnew_f_start = xnum
          if "Read Particles (New Translator) start" in line or \
              "Write Particles (New Translator) start" in line or \
              "Read Removed Particles (New Translator) start" in line or \
              "Write Removed Particles (New Translator) start" in line or \
              "Delete Particles (New Translator) start" in line or \
              "Delete Removed Particles (New Translator) start" in line:
              tnew_p_start = xnum
          if "Read Scalar Array (Old Translator) start" in line or \
              "Write Scalar Array (Old Translator) start" in line:
              told_s_start = xnum
          if "Read Face Array (Old Translator) start" in line or \
              "Write Face Array (Old Translator) start" in line:
              told_f_start = xnum
          if "Read Particles (Old Translator) start" in line or \
              "Write Particles (Old Translator) start" in line or \
              "Read Removed Particles (Old Translator) start" in line or \
              "Write Removed Particles (Old Translator) start" in line or \
              "Delete Particles (Old Translator) start" in line or \
              "Delete Removed Particles (Old Translator) start" in line:
              told_p_start = xnum
          if "Completed application" in line:
              worker_end = xnum
              worker_total += worker_end - worker_start
          if "App compute job end" in line:
              if 'loop' in line:
                  compute_non_sterile_end = xnum
                  compute_non_sterile_total += compute_non_sterile_end - compute_non_sterile_start
                  compute_non_sterile = False
              else:
                  compute_sterile_end = xnum
                  compute_sterile_total += compute_sterile_end - compute_sterile_start
                  compute_sterile = False
                  if 'projection' in line:
                      projection_end = xnum
                      projection_total += projection_end - projection_start
                      projection = False
          if "App copy job end" in line:
              copy_end = xnum
              copy_total += copy_end - copy_start
          if "Read Scalar Array (New Translator) end" in line or \
              "Write Scalar Array (New Translator) end" in line:
              tnew_s_end = xnum
              if compute_non_sterile:
                tnew_s_total_compute_non_sterile += tnew_s_end - tnew_s_start
              elif compute_sterile:
                tnew_s_total_compute_sterile += tnew_s_end - tnew_s_start
                if projection:
                  projection_trans += tnew_s_end - tnew_s_start
              elif copy:
                tnew_s_total_copy += tnew_s_end - tnew_s_start
              else:
                unknown_trans += tnew_s_end - tnew_s_start
          if "Read Face Array (New Translator) end" in line or \
              "Write Face Array (New Translator) end" in line:
              tnew_f_end = xnum
              if compute_non_sterile:
                tnew_f_total_compute_non_sterile += tnew_f_end - tnew_f_start
              elif compute_sterile:
                tnew_f_total_compute_sterile += tnew_f_end - tnew_f_start
                if projection:
                  projection_trans += tnew_f_end - tnew_f_start
              elif copy:
                tnew_f_total_copy += tnew_f_end - tnew_f_start
              else:
                unknown_trans += tnew_f_end - tnew_f_start
          if "Read Particles (New Translator) end" in line or \
              "Write Particles (New Translator) end" in line or \
              "Readgt Removed Particles (New Translator) end" in line or \
              "Write Removed Particles (New Translator) end" in line or \
              "Delete Particles (New Translator) end" in line or \
              "Delete Removed Particles (New Translator) end" in line:
              tnew_p_end = xnum
              if compute_non_sterile:
                tnew_p_total_compute_non_sterile += tnew_p_end - tnew_p_start
              elif compute_sterile:
                tnew_p_total_compute_sterile += tnew_p_end - tnew_p_start
                if projection:
                  projection_trans += tnew_p_end - tnew_p_start
              elif copy:
                tnew_p_total_copy += tnew_p_end - tnew_p_start
              else:
                unknown_trans += tnew_p_end - tnew_p_start
          if "Read Scalar Array (Old Translator) end" in line or \
              "Write Scalar Array (Old Translator) end" in line:
              told_s_end = xnum
              if compute_non_sterile:
                told_s_total_compute_non_sterile += told_s_end - told_s_start
              elif compute_sterile:
                told_s_total_compute_sterile += told_s_end - told_s_start
                if projection:
                  projection_trans += told_s_end - told_s_start
              elif copy:
                told_s_total_copy += told_s_end - told_s_start
              else:
                unknown_trans += told_s_end - told_s_start
          if "Read Face Array (Old Translator) end" in line or \
              "Write Face Array (Old Translator) end" in line:
              told_f_end = xnum
              if compute_non_sterile:
                told_f_total_compute_non_sterile += told_f_end - told_f_start
              elif compute_sterile:
                told_f_total_compute_sterile += told_f_end - told_f_start
                if projection:
                  projection_trans += told_f_end - told_f_start
              elif copy:
                told_f_total_copy += told_f_end - told_f_start
              else:
                unknown_trans += told_f_end - told_f_start
          if "Read Particles (Old Translator) end" in line or \
              "Write Particles (Old Translator) end" in line or \
              "Read Removed Particles (Old Translator) end" in line or \
              "Write Removed Particles (Old Translator) end" in line or \
              "Delete Particles (Old Translator) end" in line or \
              "Delete Removed Particles (Old Translator) end" in line:
              told_p_end = xnum
              if compute_non_sterile:
                told_p_total_compute_non_sterile += told_p_end - told_p_start
              elif compute_sterile:
                told_p_total_compute_sterile += told_p_end - told_p_start
                if projection:
                  projection_trans += told_p_end - told_p_start
              elif copy:
                told_p_total_copy += told_p_end - told_p_start
              else:
                unknown_trans += told_p_end - told_p_start
      if num == total_lines/4:
          print "Parsed 25% ..."
      if num == total_lines/2:
          print "Parsed 50% ..."
      if num == 3*total_lines/4:
          print "Parsed 75% ..."

  if worker_total == 0:
      worker_total = latest_seen - worker_start

  idle_time = worker_total - copy_total - compute_non_sterile_total - compute_sterile_total - unknown_trans

  compute_non_sterile_trans =  \
                   tnew_s_total_compute_non_sterile + tnew_f_total_compute_non_sterile + tnew_p_total_compute_non_sterile \
                 + told_s_total_compute_non_sterile + told_f_total_compute_non_sterile + told_p_total_compute_non_sterile

  compute_sterile_trans =  \
                   tnew_s_total_compute_sterile + tnew_f_total_compute_sterile + tnew_p_total_compute_sterile \
                 + told_s_total_compute_sterile + told_f_total_compute_sterile + told_p_total_compute_sterile

  copy_trans =  \
                   tnew_s_total_copy + tnew_f_total_copy + tnew_p_total_copy \
                 + told_s_total_copy + told_f_total_copy + told_p_total_copy

  compute_non_sterile_bare = compute_non_sterile_total - compute_non_sterile_trans

  compute_sterile_bare = compute_sterile_total - compute_sterile_trans

  copy_bare = copy_total - copy_trans

  projection_bare = projection_total - projection_trans

  print "Parsed 100% ..."

  
  print "Opening %s ..." % data_file
  data = open(data_file, 'w')
  data.write("Application %2.2f\n" % worker_total)
  data.write("Compute Non-Sterile Total %2.2f\n" % compute_non_sterile_total)
  data.write("Compute Sterile Total %2.2f\n" % compute_sterile_total)
  data.write("Copy Total %2.2f\n" % copy_total)

  data.write("Translator_New_Scalar From Compute Non-Sterile %2.2f\n" % tnew_s_total_compute_non_sterile)
  data.write("Translator_New_Face From Compute Non-Sterile %2.2f\n" % tnew_f_total_compute_non_sterile)
  data.write("Translator_New_Particles From Compute Non-Sterile %2.2f\n" % tnew_p_total_compute_non_sterile)
  data.write("Translator_Old_Scalar From Compute Non-Sterile %2.2f\n" % told_s_total_compute_non_sterile)
  data.write("Translator_Old_Face From Compute Non-Sterile %2.2f\n" % told_f_total_compute_non_sterile)
  data.write("Translator_Old_Particles From Compute Non-Sterile %2.2f\n" % told_p_total_compute_non_sterile)

  data.write("Translator_New_Scalar From Compute Sterile %2.2f\n" % tnew_s_total_compute_sterile)
  data.write("Translator_New_Face From Compute Sterile %2.2f\n" % tnew_f_total_compute_sterile)
  data.write("Translator_New_Particles From Compute Sterile %2.2f\n" % tnew_p_total_compute_sterile)
  data.write("Translator_Old_Scalar From Compute Sterile %2.2f\n" % told_s_total_compute_sterile)
  data.write("Translator_Old_Face From Compute Sterile %2.2f\n" % told_f_total_compute_sterile)
  data.write("Translator_Old_Particles From Compute Sterile %2.2f\n" % told_p_total_compute_sterile)

  data.write("Translator_New_Scalar From Copy %2.2f\n" % tnew_s_total_copy)
  data.write("Translator_New_Face From Copy %2.2f\n" % tnew_f_total_copy)
  data.write("Translator_New_Particles From Copy %2.2f\n" % tnew_p_total_copy)
  data.write("Translator_Old_Scalar From Copy %2.2f\n" % told_s_total_copy)
  data.write("Translator_Old_Face From Copy %2.2f\n" % told_f_total_copy)
  data.write("Translator_Old_Particles From Copy %2.2f\n" % told_p_total_copy)

  data.write("Translator From Unknown Source %2.2f\n" % unknown_trans)

  data.write("Decomposed Compute Non-Sterile Bare  %2.2f\n" % compute_non_sterile_bare)
  data.write("Decomposed Compute Non-Sterile Translator %2.2f\n" % compute_non_sterile_trans)
  data.write("Decomposed Compute Sterile Bare  %2.2f\n" % compute_sterile_bare)
  data.write("Decomposed Compute Sterile Translator %2.2f\n" % compute_sterile_trans)
  data.write("Decomposed Copy Bare  %2.2f\n" % copy_bare)
  data.write("Decomposed Copy Translator %2.2f\n" % copy_trans)
  data.write("Decomposed Idle  %2.2f\n" % idle_time)
  data.write("Decomposed Projection Bare  %2.2f\n" % projection_bare)
  data.write("Decomposed Projection Translator  %2.2f\n" % projection_trans)

  data.close()
  log.close()
  log_temp.close()


