#!/usr/bin/env python

import sys
import operator

def IsCompute(job_name):
    return job_name[0:8] == "Compute:"

def break_line(line):
    contents = line.split(' ')
    timestamp = float(contents[0])
    event = contents[1]
    if event in {"local_done", "io_done"}:
        return (None, None, None, False)
    job_name = contents[2]
    if not IsCompute(job_name):
        return (None, None, None, False)
    job_id = int(contents[3])
    name_map = {"f":"fin", "r":"run", "recv_job":"recv",
            "dispatch_job(new)":"dispatch",
            "dispatch_job(job_done)":"dispatch"}
    event = name_map[event]
    return (timestamp, event, job_id, True)


N = int(sys.argv[1])
raw_file = sys.argv[2]
result_file = sys.argv[3]
blockingtime_file = sys.argv[4]

for rank in range(1, N+1):
    print rank

    g = open(result_file.format(rank), "w")

    f = open(raw_file.format(rank), "r")
    line = f.readline()
    all_jobs = set()
    while line:
        timestamp, event, job_id, flag = break_line(line)
        if flag:
            all_jobs.add(job_id)
            g.write("{!r} {} {}\n".format(timestamp, event, job_id))
        line = f.readline()
    f.close()

    f = open(blockingtime_file, "r")
    line = f.readline()
    while line:
        contents = line.split(" ")
        job_id, compute_unblock_ts, io_unblock_ts = int(contents[0]), float(contents[1]), float(contents[2])
        if job_id in all_jobs:
            g.write("{!r} {} {}\n".format(compute_unblock_ts, "unblock_comp", job_id))
            g.write("{!r} {} {}\n".format(io_unblock_ts, "unblock_io", job_id))
        line = f.readline()
    f.close()
    g.close()
