#!/usr/bin/env python

import sys
import re

def break_line(line):
    contents = line.split(' ')
    timestamp = float(contents[0])
    event = contents[1]
    job_name = contents[2]
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
backend_log_file=sys.argv[2]
beforeset_log_file=sys.argv[3]
blockingtime_log_file=sys.argv[4]

# load
finish_dict=dict()
for rank in range(1, N+1):
    print "loading", rank
    f = open(backend_log_file.format(rank), "r")
    line = f.readline()
    while line:
        timestamp, event, job_name, job_id = break_line(line)
        line = f.readline()
        if event == 'f':
            finish_dict[job_id] = (timestamp, job_name)
    f.close()

g = open(blockingtime_log_file, "w")
f = open(beforeset_log_file, "r")
line = f.readline()
while line:
    job_id, before_set = extract_set(line)
    if job_id < 10000000000:
        line = f.readline()
        continue
    try:
        max_job_id = None
        for i_job in before_set:
            if finish_dict[i_job][1][:6] == "Comput":
                if not max_job_id or finish_dict[i_job][0] > finish_dict[max_job_id][0]:
                    max_job_id = i_job
        last_timestamp = finish_dict[max_job_id][0]
        last_name = finish_dict[max_job_id][1]
        if not job_id in finish_dict:
            raise KeyError
    except KeyError:
        pass
    else:
        # job blocked by job(id name timestamp)
        g.write("{} {} {!r} {} {} {!r}\n".format(job_id, finish_dict[job_id][1], finish_dict[job_id][0], max_job_id, last_name, last_timestamp))
    line = f.readline()
g.close()
f.close()
