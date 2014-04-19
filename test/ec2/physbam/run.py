#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import sys
import os
import time
import subprocess

import physbam_utils

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import utils
import config




utils.run_instances();
utils.wait_for_instances_to_start();

ip_addresses = utils.get_ip_addresses();
physbam_utils.copy_nodes_file_to_hosts(ip_addresses)
physbam_utils.run_experiment(ip_addresses[0])
physbam_utils.collect_output_data(ip_addresses)



utils.terminate_instances();
