#!/usr/bin/env python

import decimal
from optparse import OptionParser
import os
import os.path
import re


import numpy as np
import matplotlib.pyplot as plt
from operator import add

N = 3
t1 = (1870, 2340, 0)
p   = (1870, 2951, 0)
t2 = (199, 149, 0)
c = (344, 676, 0)
idle = (2596, 1712, 0)


t1p = map(add, t1, p)
t1pt2 = map(add, t1p, t2)
t1pt2c = map(add, t1pt2, c)


ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence

p1 = plt.bar(ind, t1, width, color='r')
p2 = plt.bar(ind, p , width, color='y', bottom=t1)
p3 = plt.bar(ind, t2, width, color='b', bottom=t1p)
p4 = plt.bar(ind, c , width, color='g', bottom=t1pt2)
p5 = plt.bar(ind, idle , width, color='w', bottom=t1pt2c)

plt.ylabel('Time (seconds)')
plt.title('Nimbus')
plt.xticks(ind+width/2., ('W1', 'W4', '') )
plt.yticks(np.arange(0,81,10))
plt.legend( (p1[0], p2[0], p3[0], p4[0], p5[0]), ('Translator Compute', 'Compute', 'Translator Copy', 'Copy', 'Idle') )

plt.show()
