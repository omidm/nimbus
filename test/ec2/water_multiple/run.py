#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import sys
import os
import time
import subprocess

import utils
import config

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ec2

# ec2.run_instances(
#     config.EC2_LOCATION,
#     config.NIMBUS_0_3_AMI,
#     config.INSTANCE_NUM,
#     config.KEY_NAME,
#     config.SECURITY_GROUP,
#     config.PLACEMENT_GROUP,
#     config.INSTANCE_TYPE);

# ec2.wait_for_instances_to_start(
#     config.EC2_LOCATION,
#     config.INSTANCE_NUM);

ip_addresses = ec2.get_ip_addresses(config.EC2_LOCATION);

scheduler_ip = "54.185.152.63"
worker_ips = list(ip_addresses)
worker_ips.remove(scheduler_ip)
# scheduler_ip = ip_addresses[0]
# worker_ips = list(ip_addresses)
# worker_ips.pop(0)

print "scheduler IP: " + scheduler_ip
print "Worker IPs: "
print worker_ips

# utils.test_workers(worker_ips)
 
# utils.build_binaries(scheduler_ip);
# utils.distribute_binaries(scheduler_ip, worker_ips);
 
# utils.run_experiment(scheduler_ip, worker_ips)
# utils.collect_output_data(scheduler_ip, worker_ips)
# utils.terminate_experiment(scheduler_ip, worker_ips)
# utils.clean_output_data(scheduler_ip, worker_ips)
 
# ec2.terminate_instances(config.EC2_LOCATION);
