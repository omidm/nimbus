#!/usr/bin/env python

import sys
import os.path
import re

def PrintErrorAndExit():
    print "Usage: ./generate_data.py <data config file>"
    sys.exit(1)

if len(sys.argv) != 2:
    print "Expected one argument, got " + str(len(sys.argv)-1)
    PrintErrorAndExit()

data_config_file = sys.argv[1]
if not os.path.isfile:
    print "Could not find file " + data_config_file
    PrintErrorAndExit()

def ValidateSizeTuples(sizes, num):
    if len(sizes) != 3:
        print "Expected 3 tuples - domain, number of partitions, ghost width"
        print "Got " + str(len(sizes)) + " at line " + str(num)
        return False
    for tup in sizes:
        if len(tup) != 3:
            print "Expect tuples of type (a, b, c)"
            print "Got " + str(tup) + " at line " + str(num)
            return False
    return True

print "Reading data configuration file " + data_config_file
data_config = open(data_config_file, 'r')

for num, line in enumerate(data_config):
    # empty lines
    if not line.split():
        continue
    # comment
    if line[0] == "#":
        continue
    # begin parsing
    args = re.split(':', line)
    if len(args) != 3:
        print "Cannot parse line " + str(num) + \
                " because it contains only " + str(len(args)-1) + " ':'"
        sys.exit(2)
    # C++ class
    cpp_class = args[0]
    # nimbus types
    nimbus_types = [x for x in re.split(' |,', args[1]) if x]
    # tuples - domain, number of partitions, ghost width
    sizes_tuple_str = [x for x in re.split('\(|\)', args[2]) if x]
    sizes_str = []
    for tup in sizes_tuple_str:
        temp = [x for x in re.split('\s|,', tup) if x]
        if (len(temp) > 0):
            sizes_str.append(temp)
    if not ValidateSizeTuples(sizes_str, num):
        sys.exit(2)
    # generate code
    print str(cpp_class) + " -- " + str(nimbus_types) + " -- " + str(sizes_str)
