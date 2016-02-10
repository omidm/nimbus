#!/usr/bin/python
import sys
import os
import time

def Loop_Watch_And_Run(command_pattern,watch_path_pattern,sequence,sleep_time):
    for i in sequence:
        done=False
        while(not done):
            watch_path=watch_path_pattern%i
            if(os.path.exists(watch_path)):
                print "path %s exists"%watch_path
                command=command_pattern%i
                print "running command",command
                os.system(command)
                done=True
            else:
                print "waiting for path %s"%watch_path
                time.sleep(sleep_time)

def main():
    if(len(sys.argv)<7):
        print "Usage: %s <command_pattern> <watch_path_pattern> <sequence_start> <sequence_stop> <sequence_step> <sleep_time>"%sys.argv[0]
        sys.exit(1)

    command_pattern=sys.argv[1]
    watch_path_pattern=sys.argv[2]
    sequence_start=int(sys.argv[3])
    sequence_stop=int(sys.argv[4])
    sequence_step=int(sys.argv[5])
    sleep_time=float(sys.argv[6])

    print "******************************************"
    print "command_pattern=",command_pattern
    print "watch_path_pattern=",watch_path_pattern
    print "sequence_start=",sequence_start
    print "sequence_stop=",sequence_stop
    print "sequence_step=",sequence_step
    print "sleep_time=",sleep_time
    print "******************************************"

    sequence=range(sequence_start,sequence_stop,sequence_step)
    Loop_Watch_And_Run(command_pattern,watch_path_pattern,sequence,sleep_time)

if __name__=="__main__":
    main()
