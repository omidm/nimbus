#!/usr/bin/env python
"""Quick summary of time parts.
"""

import sys

def is_compute(job_name):
    return len(job_name) >= 8 and job_name[0:8] == "Compute:"

def break_line(line):
    contents = line.split(' ')
    timestamp = float(contents[0])
    event = contents[1]
    if event in {"local_done", "io_done"}:
        return (timestamp, event, None, None, False)
    job_name = contents[2]
    if job_name == "JOB_DONE":
        return (timestamp, event, job_name, None, False)
    job_id = int(contents[3])
    return (timestamp, event, job_name, job_id, True)

def allocate(comp_running_jobs, non_comp_running_jobs, ready_jobs,
             blocking_jobs, cores):
    remain_cores = cores
    comp_core = min(remain_cores, comp_running_jobs)
    remain_cores -= comp_core
    non_comp_core = min(remain_cores, non_comp_running_jobs)
    remain_cores -= non_comp_core
    worker_block_core = min(remain_cores, ready_jobs)
    remain_cores -= worker_block_core
    block_core = min(remain_cores, blocking_jobs)
    remain_cores -= block_core
    idle_core = remain_cores
    return comp_core, non_comp_core, worker_block_core, block_core, idle_core

cores = int(sys.argv[1])
N = int(sys.argv[2])
file_name = sys.argv[3]
start_ts = float(sys.argv[4])
end_ts = float(sys.argv[5])

for rank in range(1, N+1):
    print 'counting for worker #{}'.format(rank)

    f = open(file_name.format(rank), 'r')
    line = f.readline()
    blocked_jobs = 0
    blocked_time = 0
    ready_jobs = 0
    ready_time = 0
    comp_running_jobs = 0
    comp_time = 0
    non_comp_running_jobs = 0
    non_comp_time = 0
    idle_time = 0
    comp_core, non_comp_core, worker_block_core, block_core, idle_core = \
            0, 0, 0, 0, 0
    timestamp, event, job_name, job_id, flag = break_line(line)
    last_timestamp = timestamp

    while line:
        timestamp, event, job_name, job_id, flag = break_line(line)
        line = f.readline()
        if end_ts and timestamp >= end_ts:
            break
        if not flag:
            continue
        if event == "recv_job":
            blocked_jobs += 1
        elif event == "dispatch_job(new)":
            blocked_jobs -= 1
            ready_jobs += 1
        elif event == "dispatch_job(job_done)":
            blocked_jobs -= 1
            ready_jobs += 1
        elif event == "r":
            ready_jobs -= 1
            if is_compute(job_name):
                comp_running_jobs += 1
            else:
                non_comp_running_jobs += 1
        elif event == "f":
            if is_compute(job_name):
                comp_running_jobs -= 1
            else:
                non_comp_running_jobs -= 1
        else:
            assert False
        if not start_ts or start_ts < last_timestamp:
            comp_time += comp_core * (timestamp - last_timestamp)
            non_comp_time += non_comp_core * (timestamp - last_timestamp)
            ready_time += worker_block_core * (timestamp - last_timestamp)
            blocked_time += block_core * (timestamp - last_timestamp)
            idle_time += idle_core * (timestamp - last_timestamp)
            comp_core, non_comp_core, worker_block_core, block_core, \
                idle_core = allocate(comp_running_jobs, non_comp_running_jobs,
                                     ready_jobs, blocked_jobs, cores)
        last_timestamp = timestamp

    print comp_time, non_comp_time, ready_time, blocked_time, idle_time

