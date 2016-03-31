#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import sys
import os
import time
import subprocess
import argparse

import ec2
import utils
import config


parser = argparse.ArgumentParser(description='Nimbus EC2 Manager.')
parser.add_argument(
    "-l", "--launch",
    dest="launch",
    action="store_true",
    help="launch ec2 instances")
parser.add_argument(
    "-t", "--terminate",
    dest="terminate",
    action="store_true",
    help="terminate ec2 instances")
parser.add_argument(
    "-m", "--monitor",
    dest="monitor",
    action="store_true",
    help="monitor ec2 instances")
parser.add_argument(
    "-s", "--start",
    dest="start",
    action="store_true",
    help="start the experiment")
parser.add_argument(
    "-e", "--end",
    dest="end",
    action="store_true",
    help="end the experiment")
parser.add_argument(
    "-d", "--download",
    dest="download",
    action="store_true",
    help="download the logs")
parser.add_argument(
    "-c", "--clean",
    dest="clean",
    action="store_true",
    help="clean the logs")
parser.add_argument(
    "-p", "--print",
    dest="printip",
    action="store_true",
    help="print controller and workers ip addresses")
parser.add_argument(
    "-w", "--wake_up",
    dest="wakeup",
    action="store_true",
    help="ssh test in to nodes for testing")

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
    "-pip", "--use_private",
    dest="useprivate",
    action="store_true",
    help="if specified will use the private ips for inter node communications")

args = parser.parse_args()


if (args.monitor):
  ec2.wait_for_instances_to_start(
      config.EC2_LOCATION,
      config.CONTROLLER_NUM + config.WORKER_NUM,
      placement_group=config.PLACEMENT_GROUP);

elif (args.launch):
  ans = raw_input("Are you sure you want to launch {} ec2 instances? (Enter 'yes' to proceed): ".format(config.CONTROLLER_NUM + config.WORKER_NUM))
  if (ans != 'yes'):
    print "Aborted"
    exit(0)

  print "Launching the instances ..."
  ec2.run_instances(
      config.EC2_LOCATION,
      config.NIMBUS_AMI,
      config.WORKER_NUM,
      config.KEY_NAME,
      config.SECURITY_GROUP,
      config.PLACEMENT,
      config.PLACEMENT_GROUP,
      config.WORKER_INSTANCE_TYPE);
  ec2.run_instances(
      config.EC2_LOCATION,
      config.NIMBUS_AMI,
      config.CONTROLLER_NUM,
      config.KEY_NAME,
      config.SECURITY_GROUP,
      config.PLACEMENT,
      config.PLACEMENT_GROUP,
      config.CONTROLLER_INSTANCE_TYPE);

elif (args.terminate):
  ans = raw_input("Are you sure you want to terminate all ec2 instances? (Enter 'yes' to proceed): ")
  if (ans != 'yes'):
    print "Aborted"
    exit(0)

  print "Terminating the instances ..."
  ec2.terminate_instances(
      config.EC2_LOCATION,
      placement_group=config.PLACEMENT_GROUP);

elif (args.start or args.end or args.download or args.clean or args.printip or args.wakeup):

  ip_addresses = ec2.get_ip_addresses(
      config.EC2_LOCATION,
      placement_group=config.PLACEMENT_GROUP);
  
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
  
  if (args.printip):
    print "Controller IP:         " + controller_ip
    print "Controller Private IP: " + controller_p_ip
    print "Worker IPs:           " + str(worker_ips)
    print "Worker Private IPs:   " + str(worker_p_ips)
  
  if (args.wakeup):
    utils.test_nodes(worker_ips + [controller_ip])
   
  if (args.start):
    utils.start_experiment(controller_ip, controller_p_ip, worker_ips, worker_p_ips)
  
  if(args.download):
    utils.collect_logs(controller_ip, worker_ips)
  
  if (args.end):
    utils.stop_experiment(controller_ip, worker_ips)

  if (args.clean):
    utils.clean_logs(controller_ip, worker_ips)
  
else :
  print "\n** Provide an action to perform!\n"
  print parser.print_help();


