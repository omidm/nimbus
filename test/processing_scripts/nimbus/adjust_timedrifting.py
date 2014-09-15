#!/usr/bin/env python

import sys
import os

N = int(sys.argv[1])
drift_file = sys.argv[2]
input_dir = sys.argv[3]
output_dir = sys.argv[4]

f = open(drift_file, "r")
for i in range(1, N+1):
	line = f.readline()
	offset = float(line)
	g = open((input_dir+"/{}_event_fe.txt").format(i), "r")
	h = open((output_dir+"/{}_event_fe.txt").format(i), "w")
	for line in g.readlines():
		space_index = line.find(" ")
		h.write("{!r} {}".format(
			float(line[:space_index]) - offset, line[space_index+1:]))
	g.close()
	h.close()
f.close()
f = open(drift_file, "r")
for i in range(1, N+1):
	line = f.readline()
	offset = float(line)
	g = open((input_dir+"/{}_event_be.txt").format(i), "r")
	h = open((output_dir+"/{}_event_be.txt").format(i), "w")
	for line in g.readlines():
		space_index = line.find(" ")
		h.write("{!r} {}".format(
			float(line[:space_index]) - offset, line[space_index+1:]))
	g.close()
	h.close()
f.close()
