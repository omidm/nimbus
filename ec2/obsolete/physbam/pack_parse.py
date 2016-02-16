import sys

pack_time = 0
unpack_time = 0

pack_bytes = 0
unpack_bytes = 0

put_time = 0
fill_time = 0

for i in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    fname = "pack_mpi" + str(i) + ".log"
    print fname
    f = open(fname)

    start_event = None
    start_time = None

    put_event = None
    put_start_time = None

    fill_event = None
    fill_start_time = None

    stack = []
    fill_stack = []
    
    num_bytes = None

    for line in f:
        if "Fill" in line and "PUT" in line:
            tokens = line.split()
            if put_event is None:
                put_event = tokens[2]
                put_start = tokens[0]
            else:
                if put_event != tokens[2]:
                    print "EXITING PUT"
                    exit()
                diff_time = float(tokens[0]) - float(put_start)
                put_time += diff_time
                put_event = None
                put_start = None
            continue
        if "Fill" in line:
            tokens = line.split()
            time = tokens[0]
            event = tokens[2]
            if fill_event == None:
                if tokens[-1] != "START":
                    exit()
                fill_start_time = time
                fill_event = event
            elif event == fill_event:
                if tokens[-1] != "END":
                    exit()
                if len(fill_stack) != 0:
                    exit()
                end_time = time
                elapsed_time = float(end_time) - float(fill_start_time)
                # print event + ": " + str(elapsed_time)
                fill_time += elapsed_time
                fill_start_time = None
                fill_event = None
            elif tokens[-1] == 'START':
                fill_stack.append(event)
            elif tokens[-1] == 'END':
                head = fill_stack.pop()
                if head != event:
                    exit()
            continue
        tokens = line.split()
        time = tokens[0]
        event = tokens[2] + tokens[4]
        if "START" in line and num_bytes is None:
            num_bytes = tokens[6]
        if start_event == None:
            if tokens[3] != "START":
                print tokens[3]
                print "EXITING"
                exit()
            start_time = time
            start_event = event
        elif event == start_event:
            if tokens[3] != "END":
                print tokens[3]
                print "EXITING 2"
                exit()
            if len(stack) != 0:
                print "EXITING 4"
                exit()
            end_time = time
            elapsed_time = float(end_time) - float(start_time)
            if "UNPACK" in event:
                unpack_time += elapsed_time
                unpack_bytes += int(num_bytes)
                # print unpack_bytes
            else:
                pack_time += elapsed_time
                pack_bytes += int(num_bytes)
            start_time = None
            start_event = None
            num_bytes = None
        elif tokens[3] == 'START':
            stack.append(event)
        elif tokens[3] == 'END':
            head = stack.pop()
            if head != event:
                print "EXITING 3"
                exit()

print "Pack Time: " + str(pack_time)
print "Unpack Time: " + str(unpack_time)
print
print "Put Time: " + str(put_time)
print
print "Fill Time: " + str(fill_time)
print
print "Pack + Put: " + str(pack_time + unpack_time + put_time)
print
print "Pack + Fill: " + str(pack_time + unpack_time + fill_time)
print
print "Pack Bytes: " + str(pack_bytes)
print "Unpack Bytes: " + str(unpack_bytes)




    
