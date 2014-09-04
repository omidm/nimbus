#!/usr/bin/env python

import sys

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


cores = int(sys.argv[1])
f = open(sys.argv[2], 'r')
g = open(sys.argv[3], 'w')
line = f.readline()

ready_jobs = 0
running_jobs = 0

timestamp, event, job_name, job_id, blocking_job_id = break_line(line)

idle = True
last_idle_timestamp = timestamp
while line:
	timestamp, event, job_name, job_id, blocking_job_id = break_line(line)
	# print timestamp, event, job_name, job_id, blocking_job_id 
	# event in [recv_job, dispatch_job(new), dispatch_job(job_done), r, f]
	if event == "recv_job":
		pass
	elif event == "dispatch_job(new)":
		ready_jobs += 1
		if idle:
			g.write("blame_scheduler %d %s %f %f\n" % 
				(job_id, job_name, timestamp, last_idle_timestamp))
		# blame job_id, it could be ran earlier at last_idle_timestamp.
		# the dependency-constraint time or the resourse constraint time.
	elif event == "dispatch_job(job_done)":
		ready_jobs += 1
		if idle:
			g.write("blame_load_balance %d %s %f %f %d\n" % 
				(job_id, job_name, timestamp, last_idle_timestamp, blocking_job_id))
		# blame blocking_job_id, job_id could be ran earlier.
		# the dependency-constraint time or the resourse constraint time.
	elif event == "r":
		running_jobs += 1
		ready_jobs -= 1
	elif event == "f":
		if job_name[0:8] == "Compute:":
			running_jobs -= 1
	else:
		assert false

	if ready_jobs + running_jobs < cores:
		if not idle:
			idle = True
			last_idle_timestamp = timestamp
	else:
		idle = False

	line = f.readline()

g.close()
