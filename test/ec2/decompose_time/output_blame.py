#!/usr/bin/env python

import sys

["Compute:loop_frame", "Compute:projection_main", "Compute:loop_iteration", "Compute:loop_iteration_part_two", "Compute:projection_loop_iteration"]

def break_line(line):
	contents = line.split(' ')
	timestamp = float(contents[0])
	assert timestamp >= break_line.last_timestamp
	break_line.last_timestamp = timestamp
	event = contents[1]
	job_name = contents[2]
	job_id = int(contents[3])
	if event == "dispatch_job(job_done)":
		blocking_job_id = int(contents[4])
	else:
		blocking_job_id = -1
	return (timestamp, event, job_name, job_id, blocking_job_id)

break_line.last_timestamp = -1

def IsCompute(job_name):
	return len(job_name) >= 8 and job_name[0:8] == "Compute:"

cores = int(sys.argv[1])
f = open(sys.argv[2], 'r')
g = open(sys.argv[3], 'w')
line = f.readline()

ready_jobs = 0
running_jobs = 0

timestamp, event, job_name, job_id, blocking_job_id = break_line(line)

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
while line:
	timestamp, event, job_name, job_id, blocking_job_id = break_line(line)
	# print timestamp, event, job_name, job_id, blocking_job_id 
	# event in [recv_job, dispatch_job(new), dispatch_job(job_done), r, f]
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
		if IsCompute(job_name):
			running_jobs -= 1
			try:
				running_set.remove(job_id)
			except KeyError:
				break
	else:
		assert false
	
	used_cores = min(cores, ready_jobs + running_jobs)
	blocked_cores = min(cores - used_cores, blocking_sum)
	time_blocking += (timestamp - last_timestamp) * blocked_cores
	time_empty += (timestamp - last_timestamp) * (cores - blocked_cores - used_cores) 
	time_active += (timestamp - last_timestamp) * used_cores 
	time_wall += timestamp - last_timestamp
	last_timestamp = timestamp
	line = f.readline()

g.write("time_blocking {} time_empty {} time_active {} time_wall {}\n".format(time_blocking, time_empty, time_active, time_wall))
g.close()
