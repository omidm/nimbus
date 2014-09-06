#!/usr/bin/env python

import sys
import re

expr = re.compile("[^0-9]+")
def extract(line):
	ls = re.split(expr, line)
	return int(ls[1])

def extract_set(line):
	ls = re.split(expr, line)
	t = set()
	for data in ls[2:]:
		if data:
			t.add(int(data))
	return t
	
def break_line(line):
	contents = line.split(' ')
	timestamp = float(contents[0])
	event = contents[1]
	job_name = contents[2]
	job_id = int(contents[3])
	return (timestamp, event, job_name, job_id)

def break_blame_line(line):
	contents = line.split(' ')
	event = contents[0]
	job_id = int(contents[1])
	job_name = contents[2]
	start_ts = float(contents[3])
	resource_constraint_ts = float(contents[4])
	if event == "blame_load_balance":
		blocking_job_id = int(contents[5])
	else:
		blocking_job_id = 0
	return (event, job_id, job_name, start_ts, resource_constraint_ts, blocking_job_id)

finish_time = dict()

# Load all the finish time.
n = int(sys.argv[1])
backend_file = sys.argv[2]
step_2_file = sys.argv[3]
before_set_file = sys.argv[4]
out_file = sys.argv[5]
for i in range(1, n+1):
	print i
	f = open(backend_file % i, 'r')
	for line in f.readlines():
		timestamp, event, job_name, job_id = break_line(line)
		if event == "f":
			finish_time[job_id] = (timestamp, i, job_name)
	f.close()

temp_g = open(before_set_file, 'r')
known_jobs = {extract(line) for line in temp_g.readlines()}
for rank in range(1, n+1):
	print rank
	f = open(step_2_file % rank, 'r')
	g = open(before_set_file, 'r')
	out = open(out_file % rank, 'w')
	for line in f.readlines():
		event, job_id, job_name, start_ts, resource_constraint_ts, blocking_job_id = break_blame_line(line)
		if not (job_id in known_jobs):
			continue
		temp = g.readline()
		while temp and extract(temp) != job_id:
			temp = g.readline()
		before_set = extract_set(temp)
		# examine whether before_set is all in the finish time log.
		if blocking_job_id != 0:
			assert blocking_job_id in finish_time
			is_remote = (finish_time[blocking_job_id][1] != rank)
			is_remote_io = (finish_time[blocking_job_id][2] == "RemoteCopyReceive")
		else:
			is_remote = False
			is_remote_io = False
		if event == "blame_scheduler":
			dependency_constraint_ts = -1
			for i in before_set:
				assert i in finish_time
				if finish_time[i][2][:8] == "Compute:":
					dependency_constraint_ts = max(dependency_constraint_ts, finish_time[i][0])
			if dependency_constraint_ts == -1:
				print job_name
				continue
			out.write("blame_scheduler (%d %s) %f\n" % (job_id, job_name, start_ts - max(resource_constraint_ts, dependency_constraint_ts)))
		elif is_remote or is_remote_io:
			before_set.remove(blocking_job_id)
			dependency_constraint_ts = -1
			for i in before_set:
				assert i in finish_time
				if finish_time[i][2][:8] == "Compute:" and finish_time[i][1] == rank:
					dependency_constraint_ts = max(dependency_constraint_ts, finish_time[i][0])
			assert dependency_constraint_ts != -1
			out.write("blame_load_balance (%d %s) %f (%d %s %d)\n" %
				(job_id, job_name, start_ts - max(resource_constraint_ts, dependency_constraint_ts),
				blocking_job_id, finish_time[blocking_job_id][2], finish_time[blocking_job_id][1]))
		else:
			before_set.remove(blocking_job_id)
			dependency_constraint_ts = -1
			for i in before_set:
				assert i in finish_time
				if finish_time[i][2][:8] == "Compute:":
					dependency_constraint_ts = max(dependency_constraint_ts, finish_time[i][0])
			if dependency_constraint_ts == -1:
				print job_name
				continue
			out.write("blame_worker (%d %s) %f (%d %s)\n" %
				(job_id, job_name, start_ts - max(resource_constraint_ts, dependency_constraint_ts),
				blocking_job_id, finish_time[blocking_job_id][2]))
	f.close()
	g.close()
	out.close()
