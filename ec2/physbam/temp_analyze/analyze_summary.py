#!/usr/bin/python
"""For mpi*.log analysis.
"""
import sys

def main(n, file_name, loop_name):
    handle_g = open(loop_name, 'r')
    loop_number = int(handle_g.readline())
    handle_g.close()
    for rank in range(0, n):
        handle_f = open(file_name.format(n), 'r')
        last_event = None
        last_timestamp = None
        running_time = 0
        blocking_time = 0
        for line in handle_f:
            timestamp, event = float(line.split()[0]), line.split()[1]
            if event == 'BLOCK':
                if last_event:
                    running_time += (timestamp - last_timestamp)
                last_event = event
                last_timestamp = timestamp
            elif event == 'RESTART':
                if last_event:
                    blocking_time += (timestamp - last_timestamp)
                last_event = event
                last_timestamp = timestamp
    print running_time/loop_number, blocking_time/loop_number

if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2], sys.argv[3])
