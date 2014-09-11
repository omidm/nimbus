#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
import operator

# ["Compute:loop_frame", "Compute:projection_main", "Compute:loop_iteration", "Compute:loop_iteration_part_two", "Compute:projection_loop_iteration"]


def IsCompute(job_name):
	return len(job_name) >= 8 and job_name[0:8] == "Compute:"

def break_line(line):
	contents = line.split(' ')
	timestamp = float(contents[0])
	event = contents[1]
	flag = False
	if event in {"local_done", "io_done"}:
		flag = True
		return (timestamp, event, None, None, None, flag)
	job_name = contents[2]
	if not IsCompute(job_name):
		flag = True	
		return (timestamp, event, job_name, None, None, flag)
	job_id = int(contents[3])
	if event == "dispatch_job(job_done)":
		blocking_job_id = int(contents[4])
	else:
		blocking_job_id = -1
	return (timestamp, event, job_name, job_id, blocking_job_id, flag)


cores = int(sys.argv[1])
N = int(sys.argv[2])
file_name = sys.argv[3]
end_ts = float(sys.argv[4])

def cal_color(running_jobs, ready_jobs, blocking_sum):
	#green_list = ["#edf8fb", "#b2e2e2", "#66c2a4", "#2ca25f", "#006d2c"]
	#red_list = ["#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"]
	#if running_jobs + ready_jobs:
	#	return green_list[min(running_jobs+ready_jobs, cores)]
	#if blocking_sum:
	#	return red_list[min(blocking_sum, cores)]
	#return "white"
	if running_jobs + ready_jobs > 0:
		return "green"
	if blocking_sum > 0:
		return "red"
	return "white"
	

for rank in range(1, N+1):
	print rank
	f = open(file_name.format(rank), 'r')
	line = f.readline()
	
	ready_jobs = 0
	running_jobs = 0
	timestamp, event, job_name, job_id, blocking_job_id, flag = break_line(line)

	idle = True
	last_timestamp = timestamp
	blocking_sum = 0
	blocking_compute = 0
	blocking_transfer = 0
	
	blocking_set = set()
	ready_set = set()
	running_set = set()
	
	time_active = 0
	time_blocking = 0
	time_empty = 0
	time_wall = 0
	
	history_color = "white"
	while line:
		timestamp, event, job_name, job_id, blocking_job_id, flag = break_line(line)
		line = f.readline()
		if end_ts and timestamp >= end_ts:
			break
		if flag:
			continue
		if event == "recv_job":
			blocking_sum += 1
			blocking_set.add(job_id)
		elif event == "dispatch_job(new)":
			blocking_sum -= 1
			try:
				blocking_set.remove(job_id)
			except KeyError:
				break
			ready_jobs += 1
			ready_set.add(job_id)
		elif event == "dispatch_job(job_done)":
			blocking_sum -= 1
			try:
				blocking_set.remove(job_id)
			except KeyError:
				break
			ready_jobs += 1
			ready_set.add(job_id)
		elif event == "r":
			ready_jobs -= 1
			try:
				ready_set.remove(job_id)
			except KeyError:
				break
			running_jobs += 1
			running_set.add(job_id)
		elif event == "f":
			running_jobs -= 1
			try:
				running_set.remove(job_id)
			except KeyError:
				break
		else:
			assert false
	
		plot_color = cal_color(running_jobs, ready_jobs, blocking_sum)
		if history_color != plot_color:
			plt.fill_between([last_timestamp, timestamp],
				[rank-1]*2, [rank]*2,
				facecolor=history_color, edgecolor=history_color) 
			last_timestamp = timestamp
			history_color = plot_color

plt.show()
