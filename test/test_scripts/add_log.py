#!/usr/bin/env python

import sys
import re
from decimal import *

file_name = sys.argv[1]

f = open(file_name, 'r')
content = f.readlines();

regexp = 'length\(s\):\s*(\d+\.\d+)\s*time\(s\):\s*(\d+\.+\d+)'
wall_time = 0;
execution_time = 0;
for line in content:
  result = re.findall(regexp, line)
  if (len(result) >= 1):
    execution_time += float(result[0][0])
    wall_time = float(result[0][1])
  else:
    print result

overhead = wall_time - execution_time

print "execution time: %.2f" % execution_time
print "wall time: %.2f" % wall_time
print "overhead: %.2f (%.2f%%)" % (overhead , overhead/wall_time*100)


