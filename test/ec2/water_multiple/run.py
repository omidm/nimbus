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

scheduler_ip = "54.190.117.78"

# worker_ips = ["54.71.30.252",
#               "54.190.255.94",
#               "54.190.76.192",
#               "54.70.165.138",
#               "54.189.39.139",
#               "54.70.68.14",
#               "54.188.214.119",
#               "54.189.244.232"]

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




# How to count number od iterations from std::out of controller? 
# either count loop_iteration_part_two jobs: 
#     cat log | grep complex | grep "projection_loop_iteration_end\." -c
#                      +
#     cat omid | grep Picked | grep " loop_iteration_part_two\."
#
# or count loop_iteration jobs (loop_frame is not templetized):
#    cat omid | grep complex | grep "loop_iteration_part_two\." -c
#                     +
#    cat omid | grep Picked | grep " loop_iteration\." -c
#

