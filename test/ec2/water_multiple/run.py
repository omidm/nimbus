#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import sys
import os
import time
import subprocess
import argparse

import utils
import config

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ec2

action_help = "action to do: "
action_help += "\nw: test the workers to see if they ar awake"
action_help += "\nr: run the experiment"
action_help += "\nc: collect the results"
action_help += "\nt: terminate and clean"
action_help += "\nk: terminate ec2 instances"

parser = argparse.ArgumentParser(description='Process log files.')
parser.add_argument(
    "-a", "--action",
    dest="action",
    default='t',
    required=True,
    help=action_help)
parser.add_argument(
    "-cip", "--controller_ip",
    dest="controllerip",
    default="X.X.X.X",
    help="controler ip")
parser.add_argument(
    "-cpip", "--controller_private_ip",
    dest="controllerprivateip",
    default="X.X.X.X",
    help="controler private ip")
parser.add_argument(
    "-rc", "--random_controller",
    dest="randomcontroller",
    action="store_true",
    help="if specified will pick a random controller")
parser.add_argument(
    "-pp", "--use_private",
    dest="useprivate",
    action="store_true",
    help="if specified will use the private ips for inter node communications")

args = parser.parse_args()

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

ip_addresses   = ec2.get_ip_addresses(config.EC2_LOCATION);

if (not args.randomcontroller):
  scheduler_ip   = args.controllerip
  scheduler_p_ip = args.controllerprivateip
else:
  scheduler_ip = ip_addresses["public"][0]
  scheduler_p_ip = ip_addresses["private"][0]


worker_ips = list(ip_addresses["public"])
worker_ips.remove(scheduler_ip)

if (not args.useprivate):
  scheduler_p_ip = scheduler_ip
  worker_p_ips = list(worker_ips)
else:
  worker_p_ips = list(ip_addresses["private"])
  worker_p_ips.remove(scheduler_p_ip)

print "scheduler IP:         " + scheduler_ip
print "scheduler Private IP: " + scheduler_p_ip
print "Worker IPs:           " + str(worker_ips)
print "Worker Private IPs:   " + str(worker_p_ips)


# utils.build_binaries(scheduler_ip);
# utils.distribute_binaries(scheduler_ip, worker_ips);
 
if args.action == 'w':
  utils.test_workers(worker_ips)
 
elif args.action == 'r':
  utils.run_experiment(scheduler_ip, scheduler_p_ip, worker_ips, worker_p_ips)

elif args.action == 'c':
  utils.collect_output_data(scheduler_ip, worker_ips)

elif args.action == 't':
  utils.terminate_experiment(scheduler_ip, worker_ips)
  utils.clean_output_data(scheduler_ip, worker_ips)
 
elif args.action == 'k':
  ec2.terminate_instances(config.EC2_LOCATION);

else :
  print "Unknown action: " + args.action



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

