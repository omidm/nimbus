#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

N = int(sys.argv[1])
input_file = sys.argv[2]

total_p_b = 0
total_p_c = 0
total_n_b = 0
total_n_c = 0
for rank in range(1, N+1):
	print rank
	f = open(input_file.format(rank), 'r')
	
	# Get to the plotting start.
	event = ""
	while event != "Compute:loop_frame\n":
		line = f.readline()
		timestamp, event = float(line.split(" ")[0]), line.split(" ")[1]
	last_timestamp = timestamp
	last_event = event

	iterations = 0
	in_projection = False
	p_c = 0
	p_b = 0
	n_c = 0
	n_b = 0
	line = f.readline()
	while line:
		timestamp, event = float(line.split(" ")[0]), line.split(" ")[1]
		line = f.readline()
		if event == "Compute:projection_inner\n":
			in_projection = True	
			iterations += 1
		elif event == "Compute:calculate_dt\n":
			in_projection = False
		if in_projection:
			if event == "RESTART\n":
				p_b += timestamp - last_timestamp
			elif event == "BLOCK\n":
				p_c += timestamp - last_timestamp
		else:
			if event == "RESTART\n":
				n_b += timestamp - last_timestamp
			elif event == "BLOCK\n":
				n_c += timestamp - last_timestamp
		last_timestamp = timestamp
		last_event = event
	print "Worker#{} iterations {} projection blocking {} calculation {} non-projection blocking {} calculation {}".format(rank, iterations, p_b, p_c, n_b, n_c)
	total_p_b += p_b
	total_p_c += p_c
	total_n_b += n_b
	total_n_c += n_c

print "________SUMMARY________"
total_all = total_p_b + total_p_c + total_n_b + total_n_c
print "average over {} iterations\n projection blocking: {} {:.1%}\n projection calculation: {} {:.1%}\n non-projection blocking: {} {:.1%}\n non-projection calculation: {} {:.1%}\n".format(
	iterations,
	total_p_b / iterations, total_p_b/total_all,
	total_p_c / iterations, total_p_c/total_all,
	total_n_b / iterations, total_n_b/total_all,
	total_n_c / iterations, total_n_c/total_all)
