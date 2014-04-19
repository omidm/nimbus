#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import boto.ec2
import time

import config

def run_instances():

  ec2 = boto.ec2.connect_to_region(config.EC2_LOCATION)

  ec2.run_instances(
      NIMBUS_AMI,
      min_count = config.INSTANCE_NUM,
      max_count = config.INSTANCE_NUM,
      key_name = config.KEY_NAME,
      security_groups = [config.SG_GROUP],
      instance_type = config.INSTANCE_TYPE)




def wait_for_instances_to_start():

  ec2 = boto.ec2.connect_to_region(config.EC2_LOCATION)

  ready_count = 0
  while (ready_count < config.INSTANCE_NUM):
    ready_count = 0;
    instances = ec2.get_only_instances()
    for inst in instances:
      if inst.state == 'running':
        ready_count += 1
    print "-> number of ready instances: " + str(ready_count) + " (out of " + str(config.INSTANCE_NUM) + ") ..."
    time.sleep(10)



def start_instances():

  ec2 = boto.ec2.connect_to_region(config.EC2_LOCATION)

  instances = ec2.get_only_instances()
  for inst in instances:
    ec2.start_instances(instance_ids=[inst.id])
 


def stop_instances():

  ec2 = boto.ec2.connect_to_region(config.EC2_LOCATION)

  instances = ec2.get_only_instances()
  for inst in instances:
    ec2.stop_instances(instance_ids=[inst.id])
    

def terminate_instances():

  ec2 = boto.ec2.connect_to_region(config.EC2_LOCATION)

  instances = ec2.get_only_instances()
  for inst in instances:
    ec2.terminate_instances(instance_ids=[inst.id])
    

