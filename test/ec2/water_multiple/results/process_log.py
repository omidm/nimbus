#!/usr/bin/env python

import decimal
import argparse
import os
import os.path
import re

## Parse the command line arguments ##
parser = argparse.ArgumentParser()
parser.add_argument("input", help="base name for input file name that follows with  -[worker number]")
parser.parse_args()




# parser = OptionParser()
# parser.add_option("-i", "--input", dest="infile",
#                   default="cache-log", type="string",
#                   help="input log file")
# parser.add_option("-o", "--output", dest="outfile",
#                   default="cache-log-data", type="string",
#                   help="output data from log file")
# parser.add_option("-d", "--dir", dest="dir",
#                   default="cache_80_111", type="string",
#                   help="directory for log/ data")
# (options, args) = parser.parse_args()
# log_file  = options.dir + "/" + options.infile
# data_file = options.dir + "/" + options.outfile
# 
# ### Initialize data to get from parser ###
# worker_total  = 0
# compute_total = 0
# copy_total    = 0
# tnew_s_total  = 0
# tnew_f_total  = 0
# tnew_p_total  = 0
# told_s_total  = 0
# told_f_total  = 0
# told_p_total  = 0
# worker_total  = 0
# compute_start = 0
# copy_start    = 0
# tnew_s_start  = 0
# tnew_f_start  = 0
# tnew_p_start  = 0
# told_s_start  = 0
# told_f_start  = 0
# told_p_start  = 0
# worker_end    = 0
# compute_end   = 0
# copy_end      = 0
# tnew_s_end    = 0
# tnew_f_end    = 0
# tnew_p_end    = 0
# told_s_end    = 0
# told_f_end    = 0
# told_p_end    = 0
# copy_trans    = 0
# 
# ### Open file and begin parsing ###
# print "Opening %s ..." % log_file
# log_temp = open(log_file, 'r')
# log = open(log_file, 'r')
# 
# total_lines = sum(1 for line in log_temp)
# print "Parsing %d lines ..." % total_lines
# 
# compute = False
# copy    = False
# for num, line in enumerate(log):
#     x =  re.findall('(\d+\.\d+$|\d+e-\d+$|\d+$)', line)
#     if len(x) > 0:
#         xnum = decimal.Decimal(x[0])
#         if "Worker starts" in line:
#             worker_start = xnum
#         if "App compute job start" in line:
#             compute_start = xnum
#             compute = True
#         if "App copy job start" in line:
#             copy_start = xnum
#             copy = True
#         if "Read Scalar Array (New Translator) start" in line or \
#             "Write Scalar Array (New Translator) start" in line:
#             tnew_s_start = xnum
#         if "Read Face Array (New Translator) start" in line or \
#             "Write Face Array (New Translator) start" in line:
#             tnew_f_start = xnum
#         if "Read Particles (New Translator) start" in line or \
#             "Write Particles (New Translator) start" in line or \
#             "Read Removed Particles (New Translator) start" in line or \
#             "Write Removed Particles (New Translator) start" in line or \
#             "Delete Particles (New Translator) start" in line or \
#             "Delete Removed Particles (New Translator) start" in line:
#             tnew_p_start = xnum
#         if "Read Scalar Array (Old Translator) start" in line or \
#             "Write Scalar Array (Old Translator) start" in line:
#             told_s_start = xnum
#         if "Read Face Array (Old Translator) start" in line or \
#             "Write Face Array (Old Translator) start" in line:
#             told_f_start = xnum
#         if "Read Particles (Old Translator) start" in line or \
#             "Write Particles (Old Translator) start" in line or \
#             "Read Removed Particles (Old Translator) start" in line or \
#             "Write Removed Particles (Old Translator) start" in line or \
#             "Delete Particles (Old Translator) start" in line or \
#             "Delete Removed Particles (Old Translator) start" in line:
#             told_p_start = xnum
#         if "Completed application" in line:
#             worker_end = xnum
#             worker_total += worker_end - worker_start
#         if "App compute job end" in line:
#             compute_end = xnum
#             compute_total += compute_end - compute_start
#             compute = False
#         if "App copy job end" in line:
#             copy_end = xnum
#             copy_total += copy_end - copy_start
#         if "Read Scalar Array (New Translator) end" in line or \
#             "Write Scalar Array (New Translator) end" in line:
#             tnew_s_end = xnum
#             tnew_s_total += tnew_s_end - tnew_s_start
#             if not compute:
#                 copy_trans += tnew_s_end - tnew_s_start
#         if "Read Face Array (New Translator) end" in line or \
#             "Write Face Array (New Translator) end" in line:
#             tnew_f_end = xnum
#             tnew_f_total += tnew_f_end - tnew_f_start
#             if not compute:
#                 copy_trans += tnew_f_end - tnew_f_start
#         if "Read Particles (New Translator) end" in line or \
#             "Write Particles (New Translator) end" in line or \
#             "Read Removed Particles (New Translator) end" in line or \
#             "Write Removed Particles (New Translator) end" in line or \
#             "Delete Particles (New Translator) end" in line or \
#             "Delete Removed Particles (New Translator) end" in line:
#             tnew_p_end = xnum
#             tnew_p_total += tnew_p_end - tnew_p_start
#             if not compute:
#                 copy_trans += tnew_p_end - tnew_p_start
#         if "Read Scalar Array (Old Translator) end" in line or \
#             "Write Scalar Array (Old Translator) end" in line:
#             told_s_end = xnum
#             told_s_total += told_s_end - told_s_start
#             if not compute:
#                 copy_trans += told_s_end - told_s_start
#         if "Read Face Array (Old Translator) end" in line or \
#             "Write Face Array (Old Translator) end" in line:
#             told_f_end = xnum
#             told_f_total += told_f_end - told_f_start
#             if not compute:
#                 copy_trans += told_f_end - told_f_start
#         if "Read Particles (Old Translator) end" in line or \
#             "Write Particles (Old Translator) end" in line or \
#             "Read Removed Particles (Old Translator) end" in line or \
#             "Write Removed Particles (Old Translator) end" in line or \
#             "Delete Particles (Old Translator) end" in line or \
#             "Delete Removed Particles (Old Translator) end" in line:
#             told_p_end = xnum
#             told_p_total += told_p_end - told_p_start
#             if not compute:
#                 copy_trans += told_p_end - told_p_start
#     if num == total_lines/4:
#         print "Parsed 25% ..."
#     if num == total_lines/2:
#         print "Parsed 50% ..."
#     if num == 3*total_lines/4:
#         print "Parsed 75% ..."
# print "Parsed 100% ..."
# 
# print "Opening %s ..." % data_file
# data = open(data_file, 'w')
# data.write("Application %d\n" % worker_total)
# data.write("Computation %d\n" % compute_total)
# data.write("Copy %d\n" % copy_total)
# data.write("Translator_New_Scalar %d\n" % tnew_s_total)
# data.write("Translator_New_Face %d\n" % tnew_f_total)
# data.write("Translator_New_Particles %d\n" % tnew_p_total)
# data.write("Translator_Old_Scalar %d\n" % told_s_total)
# data.write("Translator_Old_Face %d\n" % told_f_total)
# data.write("Translator_Old_Particles %d\n" % told_p_total)
# data.write("Translator_From_Copy %d\n" % copy_trans)
