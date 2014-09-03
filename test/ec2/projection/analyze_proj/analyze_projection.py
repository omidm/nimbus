#!/usr/bin/env python
import sys
import re

def AddSum(sum_s, key, value):
	if key in sum_s:
		sum_s[key] += value
	else:
		sum_s[key] = value

file_name = sys.argv[1]
f = open(file_name, 'r')
sum_dict = dict()
for line in f.readlines():
	job_name = (re.split("[\s,]", line))[2]
	time_length = float((re.findall("[0-9]+[0-9\.]*", line))[0])
	if job_name[:7] == "Compute":
		compute_job_name = job_name[7:]
	else:
		compute_job_name = ""
	AddSum(sum_dict, job_name, time_length)
for i in sum_dict:
	print i, sum_dict[i]
