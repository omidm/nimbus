#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

def break_line(line):
    contents = line.split(" ")
    return float(contents[0]), int(contents[1]), int(contents[2]), \
        int(contents[3]), int(contents[4]), int(contents[5]), int(contents[6])

def cal_color(e0, e1, e2, e3, e4, e5):
    if e1 > 0:
        return "green"
    if e0+e2 > 0:
        return "blue"
    if e5 > 0:
        return "#dd1c77"
    if e3 > 0:
        return "red"
    if e4 > 0:
        return "yellow"
    return "white"


cores = int(sys.argv[1])
N = int(sys.argv[2])
prep_file = sys.argv[3]
start_ts = float(sys.argv[4])
end_ts = float(sys.argv[5])

for rank in range(1, N+1):
    print rank
    f = open(prep_file.format(rank), 'r')
    line = f.readline()
    
    last_timestamp, _, _, _, _, _, _ = break_line(line)
    last_timestamp -= 1e-6
    history_color = "white"
    while line:
        timestamp, e1, e2, e3, e4, e5, e6 = break_line(line)
        line = f.readline()
        if end_ts and timestamp >= end_ts:
            break
        plot_color = cal_color(e1, e2, e3, e4, e5, e6)
        if history_color != plot_color:
            assert last_timestamp != timestamp
            plt.fill_between(
                [last_timestamp - start_ts if start_ts else last_timestamp,
                timestamp - start_ts if start_ts else timestamp],
                [rank-1]*2, [rank]*2,
                facecolor=history_color, edgecolor=history_color)
            last_timestamp = timestamp
            history_color = plot_color
plt.title("Nimbus running on 8 workers, each with 8 threads, scale 256, 64 app partitions, first 3 iterations.")
plt.xlabel("time(s)")
plt.ylabel("worker id")
plt.show()
