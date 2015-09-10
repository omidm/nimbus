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
action_help += "\np: just print the ip addresses"
action_help += "\nw: test the workers to see if they ar awake"
action_help += "\nr: run the experiment"
action_help += "\nc: collect the results"
action_help += "\nt: terminate and clean"
action_help += "\ns: start ec2 instances"
action_help += "\nm: monitor ec2 instances"
action_help += "\nk: terminate ec2 instances"

parser = argparse.ArgumentParser(description='Process log files.')
parser.add_argument(
    "-a", "--action",
    dest="action",
    default='p',
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
    "-pp", "--use_private",
    dest="useprivate",
    action="store_true",
    help="if specified will use the private ips for inter node communications")

args = parser.parse_args()


if args.action == 's':
  ec2.run_instances(
      config.EC2_LOCATION,
      config.NIMBUS_AMI,
      config.WORKER_INSTANCE_NUM,
      config.KEY_NAME,
      config.SECURITY_GROUP,
      config.PLACEMENT,
      config.PLACEMENT_GROUP,
      config.WORKER_INSTANCE_TYPE);
  ec2.run_instances(
      config.EC2_LOCATION,
      config.NIMBUS_AMI,
      config.CONTROLLER_INSTANCE_NUM,
      config.KEY_NAME,
      config.SECURITY_GROUP,
      config.PLACEMENT,
      config.PLACEMENT_GROUP,
      config.CONTROLLER_INSTANCE_TYPE);

elif args.action == 'k':
  ec2.terminate_instances(
      config.EC2_LOCATION,
      placement_group=config.PLACEMENT_GROUP);

elif args.action == 'm':
  ec2.wait_for_instances_to_start(
      config.EC2_LOCATION,
      config.WORKER_INSTANCE_NUM + config.CONTROLLER_INSTANCE_NUM,
      placement_group=config.PLACEMENT_GROUP);

else:

  ip_addresses = ec2.get_ip_addresses(
      config.EC2_LOCATION,
      placement_group=config.PLACEMENT_GROUP);
  
  if args.action == 'p':
    print ip_addresses
    exit(0)

  if (not args.controllerip == "X.X.X.X"):
    controller_ip   = args.controllerip
    controller_p_ip = args.controllerprivateip
  else:
    cips = ec2.get_ip_addresses(
        config.EC2_LOCATION,
        placement_group=config.PLACEMENT_GROUP,
        instance_type=config.CONTROLLER_INSTANCE_TYPE);
    controller_ip   = cips["public"][0]
    controller_p_ip = cips["private"][0]
  
  worker_ips = list(ip_addresses["public"])
  worker_ips.remove(controller_ip)
  
  if (not args.useprivate):
    controller_p_ip = controller_ip
    worker_p_ips = list(worker_ips)
  else:
    assert(not controller_p_ip == "X.X.X.X") 
    worker_p_ips = list(ip_addresses["private"])
    worker_p_ips.remove(controller_p_ip)
  
  print "Controller IP:         " + controller_ip
  print "Controller Private IP: " + controller_p_ip
  print "Worker IPs:           " + str(worker_ips)
  print "Worker Private IPs:   " + str(worker_p_ips)
  
  
  if args.action == 'w':
    utils.test_nodes(worker_ips + [controller_ip])
   
  elif args.action == 'r':
    utils.run_experiment(controller_ip, controller_p_ip, worker_ips, worker_p_ips)
  
  elif args.action == 'c':
    utils.collect_output_data(controller_ip, worker_ips)
  
  elif args.action == 't':
    utils.terminate_experiment(controller_ip, worker_ips)
    utils.clean_output_data(controller_ip, worker_ips)
  
  else :
    print "Unknown action: " + args.action





# utils.build_binaries(controller_ip);
# utils.distribute_binaries(controller_ip, worker_ips);
 


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

