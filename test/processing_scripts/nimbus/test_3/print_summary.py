#!/usr/bin/env python

import sys
def break_line(line):
    contents = line.split(" ")
    return float(contents[0]), int(contents[1]), int(contents[2]), int(contents[3]), int(contents[4]), int(contents[5]), int(contents[6])


def decompose_time(e0, e1, e2, e3, e4, e5):
    remain_cores = cores
    comp_cores = min(remain_cores, e1)
    remain_cores -= comp_cores
    copy_comp_cores = min(remain_cores, e2)
    remain_cores -= copy_comp_cores
    worker_block_cores = min(remain_cores, e0)
    remain_cores -= worker_block_cores
    sched_block_cores = min(remain_cores, e5)
    remain_cores -= sched_block_cores
    io_block_cores = min(remain_cores, e4)
    remain_cores -= io_block_cores
    comp_block_cores = min(remain_cores, e3)
    remain_cores -= comp_block_cores
    empty_cores = remain_cores
    return copy_comp_cores, comp_cores, worker_block_cores, sched_block_cores, io_block_cores, comp_block_cores, empty_cores


cores = int(sys.argv[1])
N = int(sys.argv[2])
prep_file = sys.argv[3]

s0 = 0
s1 = 0
s2 = 0
s3 = 0
s4 = 0
s5 = 0
s6 = 0
for rank in range(1, N+1):
    print rank
    f = open(prep_file.format(rank), 'r')
    line = f.readline()

    last_timestamp, _, _, _, _, _, _ = break_line(line)
    last_timestamp -= 1e-6
    history_color = "white"
    ts0, ts1, ts2, ts3, ts4, ts5, ts6 = decompose_time(0, 0, 0, 0, 0, 0)
    while line:
        timestamp, e0, e1, e2, e3, e4, e5 = break_line(line)
        line = f.readline()
        ns0, ns1, ns2, ns3, ns4, ns5, ns6 = decompose_time(e0, e1, e2, e3, e4, e5)
        s0 += ts0 * (timestamp- last_timestamp)
        s1 += ts1 * (timestamp- last_timestamp)
        s2 += ts2 * (timestamp- last_timestamp)
        s3 += ts3 * (timestamp- last_timestamp)
        s4 += ts4 * (timestamp- last_timestamp)
        s5 += ts5 * (timestamp- last_timestamp)
        s6 += ts6 * (timestamp- last_timestamp)
        ts0, ts1, ts2, ts3, ts4, ts5, st6 = ns0, ns1, ns2, ns3, ns4, ns5, ns6
        last_timestamp = timestamp
print s0, s1, s2, s3, s4, s5, s6
print "copy compute cycles {} compute cycles {} worker_block_cores {} scheduler \
blocked cycles {} io blocked cycles {} compute blocked cycles {} empty cycles \
{}".format(s0, s1, s2, s3, s4, s5, s6)
