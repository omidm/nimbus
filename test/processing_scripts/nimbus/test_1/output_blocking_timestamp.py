#!/usr/bin/env python

import sys
import re

def break_line(line):
    """
    Breaks the line of event_be.txt.
    """
    contents = line.split(' ')
    timestamp = float(contents[0])
    event = contents[1]
    job_name = contents[2]
    job_id = int(contents[3])
    return (timestamp, event, job_name, job_id)

expr = re.compile("[^0-9]+")
def extract_set(line):
    """
    Breaks the line of log_before_set. Returns the job id and its before set.
    """
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
    print "Load into memory the job finishing time of worker #{}".format(rank)
    f = open(backend_log_file.format(rank), "r")
    line = f.readline()
    while line:
        timestamp, event, job_name, job_id = break_line(line)
        line = f.readline()
        if event == 'f':
            finish_dict[job_id] = (timestamp, job_name[:6]=="Remote")
    f.close()

g = open(blockingtime_log_file, "w")
f = open(beforeset_log_file, "r")
line = f.readline()
print "Calculate the unblocking time of each computation job."
while line:
    job_id, before_set = extract_set(line)
    try:
        temp = [finish_dict[i_job][0] for i_job in before_set if not finish_dict[i_job][1]]
        last_compute_timestamp = max(temp) if temp else 0
        temp = [finish_dict[i_job][0] for i_job in before_set]
        last_all_timestamp = max(temp) if temp else 0
    except KeyError:
        pass
    else:
        g.write("{} {!r} {!r}\n".format(job_id, last_compute_timestamp, last_all_timestamp))
    line = f.readline()
g.close()
f.close()
