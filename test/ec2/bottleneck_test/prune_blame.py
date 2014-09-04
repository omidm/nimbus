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
	f = open(backend_file % i, 'r')
	for line in f.readlines():
		timestamp, event, job_name, job_id = break_line(line)
		if event == "f":
			finish_time[job_id] = (timestamp, i, job_name)
	f.close()

for rank in range(1, n+1):
	f = open(step_2_file % rank, 'r')
	g = open(before_set_file, 'r')
	out = open(out_file % rank, 'w')
	for line in f.readlines():
		event, job_id, job_name, start_ts, resource_constraint_ts, blocking_job_id = break_blame_line(line)
		temp = g.readline()
		while temp and extract(temp) != job_id:
			temp = g.readline()
		before_set = extract_set(temp)
		if event != "blame_scheduler":
			assert blocking_job_id in before_set
		before_set.discard(blocking_job_id)
		if before_set:
			dependency_constraint_ts = max([finish_time[i][0] for i in before_set])
		else:
			dependency_constraint_ts = 0
		if event == "blame_scheduler":
			out.write("blame_scheduler (%d %s) %f\n" % (job_id, job_name, start_ts - max(resource_constraint_ts, dependency_constraint_ts)))
		else:
			out.write("blame_load_balance (%d %s) %f (%d %s %s)\n" %
				(job_id, job_name, start_ts - max(resource_constraint_ts, dependency_constraint_ts),
				blocking_job_id, finish_time[blocking_job_id][2], "local" if (finish_time[blocking_job_id][1] == rank and finish_time[blocking_job_id][2] != "RemoteCopyReceive") else "remote"))
	f.close()
	g.close()
	out.close()
