#!/usr/bin/env python

import sys
import re

file_name = sys.argv[1]

f = open(file_name, 'r')
content = f.read()

pattern = re.compile('ERROR.*')
error_num = len(re.findall(pattern, content))

pattern = re.compile('WARNING.*')
warn_num = len(re.findall(pattern, content))

print (error_num + warn_num)




