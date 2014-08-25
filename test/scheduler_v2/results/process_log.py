#!/usr/bin/env python

import sys
import argparse
import os
import re
import numpy
import decimal


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
parser.add_argument(
    "-t", "--tag",
    dest="tag",
    default="assign",
    help="tag to sum up the time for")


args = parser.parse_args()

file_name = args.idir + '/' + args.ifname
print 'Opening the file ' + file_name

f = open(file_name, 'r')
content = f.readlines()

time = 0;

for line in content:
  result = re.findall('.*' + args.tag + ': (\d+(?: |$)|\d+\.\d+(?: |$)|\d+e-\d+(?: |$)|\d+\.\d+e-\d+(?: |$)).*', line)
  if len(result) > 0:
    time += decimal.Decimal(result[0])
    print result[0]

f.close()

print "Total " + args.tag + " time: " + str(time)


