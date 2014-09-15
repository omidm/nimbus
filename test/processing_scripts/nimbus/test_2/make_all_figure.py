#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

f = open("temp", "r")
mydict = eval(f.readline())
print "finish loading"
job_id = 10000664772
if "dispatch_job(new)" in mydict[job_id]:
	last_timestamp = mydict[job_id]["dispatch_job(new)"]
else:
	last_timestamp = mydict[job_id]["dispatch_job(job_done)"]
last_worker = mydict[job_id]["worker"]

calculate_t = 0
resource_t = 0
scheduler_t = 0
control_t = 0
const = 3000
iteration = -1

while job_id in mydict:
	print mydict[job_id]
	plt.fill_between(
		[mydict[job_id]["f"]+const, last_timestamp+const],
		[last_worker-1]*2, [last_worker]*2,
		facecolor="red", edgecolor="red") 
	scheduler_t += last_timestamp - mydict[job_id]["f"]
	last_timestamp = mydict[job_id]["f"]

	plt.fill_between(
		[mydict[job_id]["r"]+const, last_timestamp+const],
		[mydict[job_id]["worker"]-1]*2, [mydict[job_id]["worker"]]*2,
		facecolor="green", edgecolor="green") 
	if "calculate_dt" in mydict[job_id]["name"]:
		iteration += 1
		if iteration == 3:
			break
	if ("loop" in mydict[job_id]["name"]) or ("projection_main" in mydict[job_id]["name"]):
		control_t += last_timestamp - mydict[job_id]["r"]
	else:
		calculate_t += last_timestamp - mydict[job_id]["r"]
	if (mydict[job_id]["name"][:23]!="Compute:projection_loop" and
		mydict[job_id]["name"][:23]!="Compute:projection_step" and
		mydict[job_id]["name"][:23]!="Compute:projection_redu"):
		plt.text(mydict[job_id]["r"], mydict[job_id]["worker"]-1, mydict[job_id]["name"], rotation="vertical", fontsize=8)
	last_timestamp = mydict[job_id]["r"]
	if "dispatch_job(new)" in mydict[job_id]:
		plt.fill_between(
			[mydict[job_id]["dispatch_job(new)"]+const, last_timestamp+const],
			[mydict[job_id]["worker"]-1]*2, [mydict[job_id]["worker"]]*2,
			facecolor="blue", edgecolor="blue") 
		resource_t += last_timestamp - mydict[job_id]["dispatch_job(new)"]
		last_timestamp = mydict[job_id]["dispatch_job(new)"]
	else:
		plt.fill_between(
			[mydict[job_id]["dispatch_job(job_done)"]+const, last_timestamp+const],
			[mydict[job_id]["worker"]-1]*2, [mydict[job_id]["worker"]]*2,
			facecolor="blue", edgecolor="blue") 
		resource_t += last_timestamp - mydict[job_id]["dispatch_job(job_done)"]
		last_timestamp = mydict[job_id]["dispatch_job(job_done)"]

	if mydict[job_id]["name"] == "Compute:main":
		break
	last_worker = mydict[job_id]["worker"]
	job_id = mydict[job_id]["block_id"]
f.close()
print "normal calculation time {} control calculation time {} resource constraint blocking {} scheduler delay/io delay {}".format(calculate_t, control_t, resource_t, scheduler_t)
plt.savefig("first_three.png")
