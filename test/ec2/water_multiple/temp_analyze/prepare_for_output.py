#!/usr/bin/python
"""
Summary the outputs.
"""

for scale in [256, 344, 400, 440, 480, 512]:
    loop_file = open('../output-sep24-scale{}-experiment1/loop'.format(scale),
                     'r')
    line = loop_file.readline()
    loop_number = int(line.strip())
    loop_file.close()
    summary_file = open('summary_{}.txt'.format(scale), 'r')
    lines = summary_file.readlines()
    total = [0] * 5
    for index in range(1, 16, 2):
        for i, value in enumerate(lines[index].split()):
            total[i] += float(value)
    for i in range(0, len(total)):
        total[i] /= loop_number * 64
    print total
    print loop_number

