#!/usr/bin/env python

import sys
import re

expr = re.compile("[^0-9\.\-]+")

def get_key(line):
	ls = re.split(expr, line)
	if ls[0]:
		return float(ls[0])
	else:
		return float(ls[1])
	
fin = sys.argv[1]
fout = sys.argv[2]
f = open(fin, 'r')
lines = f.readlines()
sorted_lines = sorted(lines, key = get_key)
g = open(fout, 'w')
for line in sorted_lines:
	g.write(line)
