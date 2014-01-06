#!/usr/bin/env python

import sys
import re

file_name = sys.argv[1]

f = open(file_name, 'r')
content = f.read()

pattern = re.compile('OUTPUT: (.*)')
result = re.findall(pattern, content)

# correct_result = "-111985389, 217213596, -305294628, 362485704, "
correct_result = sys.argv[2]; 

print (str(result[0]) == correct_result)




