#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
import operator

def break_line(line):
    contents = line.split(' ')
    timestamp = float(contents[0])
    event = contents[1]
    job_id = int(contents[2])
    return (timestamp, event, job_id)

N = int(sys.argv[1])
sort_state_file = sys.argv[2]
prepare_result_file = sys.argv[3]

for rank in range(1, N+1):
    print "Preparing the data for drawing of worker#{}".format(rank)
    f = open(sort_state_file.format(rank), 'r')
    g = open(prepare_result_file.format(rank), 'w')

    line = f.readline()
    timestamp, _, _ = break_line(line)
    last_timestamp = timestamp - 1e-6
    last_record = (0, 0, 0, 0, 0, 0)

    ready_jobs = 0
    compute_running_jobs = 0
    non_compute_running_jobs = 0
    ready_set = set()
    running_set = set()
    blocking_set = dict()
    blocking_log = dict()

    while line:
        timestamp, event, job_id = break_line(line)
        line = f.readline()
        if event == "recv":
            if job_id < 10000000000:
                continue
            blocking_set[job_id] = 1
            if job_id in blocking_log:
                blocking_set[job_id] = blocking_log[job_id]
        elif event == "dispatch":
            if job_id < 10000000000:
                continue
            try:
                del blocking_set[job_id]
            except KeyError:
                break
            ready_jobs += 1
            ready_set.add(job_id)
        elif event == "run":
            if job_id < 10000000000:
                non_compute_running_jobs += 1
            else:
                ready_jobs -= 1
                try:
                    ready_set.remove(job_id)
                except KeyError:
                    break
                compute_running_jobs += 1
                running_set.add(job_id)
        elif event == "fin":
            # delete at the latest time.
            if job_id in blocking_log:
                del blocking_log[job_id]
            if job_id < 10000000000:
                non_compute_running_jobs -= 1
            else:
                compute_running_jobs -= 1
                try:
                    running_set.remove(job_id)
                except KeyError:
                    break
        elif event == "unblock_io":
            if job_id in blocking_set:
                blocking_set[job_id] = max(3, blocking_set[job_id])
            else:
                if job_id in blocking_log:
                    blocking_log[job_id] = max(3, blocking_log[job_id])
                else:
                    blocking_log[job_id] = 3
        elif event == "unblock_comp":
            if job_id in blocking_set:
                blocking_set[job_id] = max(2, blocking_set[job_id])
            else:
                if job_id in blocking_log:
                    blocking_log[job_id] = max(2, blocking_log[job_id])
                else:
                    blocking_log[job_id] = 2
        else:
            assert false
        if last_timestamp > 0 and timestamp != last_timestamp:
            g.write(repr(last_timestamp))
            g.write(" {} {} {} {} {} {}\n".format(*last_record))
        last_timestamp = timestamp
        last_record = (ready_jobs, compute_running_jobs,
                 non_compute_running_jobs,
                 sum([1 for i in blocking_set if blocking_set[i] == 1]),
                 sum([1 for i in blocking_set if blocking_set[i] == 2]),
                 sum([1 for i in blocking_set if blocking_set[i] == 3]))
    g.close()
