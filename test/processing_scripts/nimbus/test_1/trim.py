#!/usr/bin/env python

import sys
import re

f = open(sys.argv[1], 'r')
g = open(sys.argv[2], 'w')
contents = f.readlines()
for line in contents[:-1]:
	g.write(line)
g.close()
