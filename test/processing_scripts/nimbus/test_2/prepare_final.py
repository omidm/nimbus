#!/usr/bin/env python

import sys
import re

def break_blocking_log(line):
	contents = line.split(' ')
	job_id = int(contents[0])
	name = contents[1]
	timestamp = float(contents[2])
	block_job_id = int(contents[3])
	return job_id, name, timestamp, block_job_id

def break_line(line):
	contents = line.split(' ')
	timestamp = float(contents[0])
	event = contents[1]
	job_name = contents[2]
	job_id = int(contents[3])
	return (timestamp, event, job_name, job_id)

def break_fe_line(line):
	contents = line.split(' ')
	timestamp = float(contents[0])
	event = contents[1]
	if event[:12] != "dispatch_job":
		return None
	job_name = contents[2]
	if job_name[:7] != "Compute":
		return None
	job_id = int(contents[3])
	return (timestamp, event, job_name, job_id)

expr = re.compile("[^0-9]+")
def extract_set(line):
	contents = re.split(expr, line)
	t_set = set()
	for data in contents[2:]:
		if data:
			t_set.add(int(data))
	return int(contents[1]), t_set

N=int(sys.argv[1])
frontend_log_file=sys.argv[2]
backend_log_file=sys.argv[3]
blockingtime_log_file=sys.argv[4]

# load
all_dict=dict()
for rank in range(1, N+1):
	print "loading", rank
	f = open(backend_log_file.format(rank), "r")
	line = f.readline()
	while line:
		timestamp, event, job_name, job_id = break_line(line)
		line = f.readline()
		if job_name[:7] != "Compute":
			continue
		if not job_id in all_dict:
			all_dict[job_id] = dict()
		all_dict[job_id][event] = timestamp
		all_dict[job_id]["name"] = job_name
		all_dict[job_id]["worker"] = rank
	f.close()

	f = open(frontend_log_file.format(rank), "r")
	line = f.readline()
	while line:
		result = break_fe_line(line)
		line = f.readline()
		if not result: 
			continue
		timestamp, event, job_name, job_id = result
		if not job_id in all_dict:
			all_dict[job_id] = dict()
		all_dict[job_id][event] = timestamp
		all_dict[job_id]["name"] = job_name
	f.close()

f = open(blockingtime_log_file, "r")
line = f.readline()

while line:
	job_id, name, timestamp, blocking_job_id = break_blocking_log(line)
	line = f.readline()
	if not job_id in all_dict:
		all_dict[job_id] = dict()
	all_dict[job_id]["block_id"] = blocking_job_id
	all_dict[job_id]["name"] = name

f.close()

f = open("temp", "w")
f.write(repr(all_dict))
#eval
f.close()
