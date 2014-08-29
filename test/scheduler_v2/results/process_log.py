#!/usr/bin/env python

import sys
import argparse
import os
import re
import numpy
from decimal import *


## Parse the command line arguments ##
parser = argparse.ArgumentParser(description='Process log files.')
parser.add_argument(
    "-i", "--input",
    dest="ifname",
    default="log",
    help="input file name")
parser.add_argument(
    "-d", "--directory",
    dest="idir",
    default=".",
    help="directory to find the input file")
parser.add_argument(
    "-od", "--outdirectory",
    dest="odir",
    default=".",
    help="directory to dump the output file")
parser.add_argument(
    "-o", "--output",
    dest="ofname",
    default="data",
    help="output file name")

args = parser.parse_args()

file_name = args.idir + '/' + args.ifname
print 'Opening the file ' + file_name

f = open(file_name, 'r')
content = f.readlines()

versioning_time = 0;

for line in content:
  result = re.findall('.*versioning: (\d+\.\d+).*', line)
  if len(result) > 0:
    versioning_time += float(result[0])

f.close()

print "Total versioning time: " + str(versioning_time)


