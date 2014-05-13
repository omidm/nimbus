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



ec2.run_instances(
    config.EC2_LOCATION,
    config.NIMBUS_AMI,
    config.INSTANCE_NUM,
    config.KEY_NAME,
    config.SECURITY_GROUP,
    config.PLACEMENT_GROUP,
    config.INSTANCE_TYPE);

ec2.wait_for_instances_to_start(
    config.EC2_LOCATION,
    config.INSTANCE_NUM);

ip_addresses = ec2.get_ip_addresses(config.EC2_LOCATION);
print ip_addresses

utils.copy_binary_file_to_hosts(ip_addresses)
utils.copy_nodes_file_to_hosts(ip_addresses)
utils.run_experiment(ip_addresses[0])
utils.collect_output_data(ip_addresses)

ec2.terminate_instances(config.EC2_LOCATION); 
