#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

N = int(sys.argv[1])
max_iteration = int(sys.argv[2])
input_file = sys.argv[3]

for rank in range(1, N+1):
	print rank
	f = open(input_file.format(rank), 'r')
	
	# Get to the plotting start.
	event = ""
	while event != "Compute:loop_frame\n":
		line = f.readline()
		timestamp, event = float(line.split(" ")[0]), line.split(" ")[1]
	last_timestamp = timestamp

	iterations = -1
	line = f.readline()
	while line:
		timestamp, event = float(line.split(" ")[0]), line.split(" ")[1]
		if iterations == max_iteration:
			break
		line = f.readline()
		if event == "Compute:calculate_dt\n":
			iterations += 1
		plot_color = "red" if event == "RESTART\n" else "green"
		plt.fill_between(
			[last_timestamp, timestamp],
			[rank-1]*2, [rank]*2,
			facecolor=plot_color, edgecolor=plot_color) 
		last_timestamp = timestamp
plt.title("MPI on 8 machines(c3.2xlarge) and 64 processors. Scale=256. 55 iterations, 3 frames.")
plt.xlabel("time(s)")
plt.ylabel("mpi processor id")
plt.show()
