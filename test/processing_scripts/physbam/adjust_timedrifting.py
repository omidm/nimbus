#!/usr/bin/env python

import sys
import os

# M machines
M = int(sys.argv[1])
# N cores each machine
N = int(sys.argv[2])
input_file = sys.argv[3]
output_file = sys.argv[4]

constant = 100

for i in range(1, M+1):
	timestamps = []
	for j in range(N*(i-1)+1,N*i+1):
		f = open(input_file.format(j), "r")
		lines = f.readlines()
		# the line that calculate_dt is synced.
		timestamps.append(float(lines[12].split(" ")[0]))
		f.close()
	offset = sum(timestamps) / len(timestamps)

	for j in range(N*(i-1)+1,N*i+1):
		print j, repr(offset)
		f = open(input_file.format(j), "r")
		g = open(output_file.format(j), "w")
		for line in f.readlines():
			space_index = line.find(" ")
			g.write("{!r} {}".format(
				float(line[:space_index]) - offset + constant, line[space_index+1:]))
		f.close()
		g.close()
        
