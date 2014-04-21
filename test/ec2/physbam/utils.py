#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import sys
import os
import subprocess

import config

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ec2

temp_file_name = '_temp_file_'

def make_nodes_file_content(ip_addresses):

  string = ""
  for ip in ip_addresses:
    print ip
    string = string + ip + " cpu=1\n"

  file = open(temp_file_name, 'w+')
  file.write(string)
  file.close()



def copy_nodes_file_to_hosts(ip_addresses):

  make_nodes_file_content(ip_addresses)

  for ip in ip_addresses:
    subprocess.call(['scp', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        temp_file_name,
        'ubuntu@' + ip + ':~/' + config.DIRECTORY_PATH + config.NODES_FILE_NAME])

  subprocess.call(['rm', temp_file_name])


def run_experiment(ip):
  command =  'cd ' + config.DIRECTORY_PATH + ';'
  command += 'mpirun -hostfile ' + config.NODES_FILE_NAME
  command += ' -np ' + str(config.INSTANCE_NUM)
  command += ' ./Water -scale ' + str(config.SCALE)
  command += ' -e ' + str(config.FRAME_NUM)

  subprocess.call(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + ip, command])


def collect_output_data(ip_addresses):

  subprocess.call(['rm', '-rf', config.OUTPUT_NAME])
  subprocess.call(['mkdir', '-p', config.OUTPUT_NAME])

  process_num = 0
  for ip in ip_addresses:
    process_num += 1
    subprocess.call(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':~/' + config.DIRECTORY_PATH + config.OUTPUT_NAME + str(process_num),
        config.OUTPUT_NAME])



