#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import sys
import os
import subprocess

import config

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ec2

temp_file_name = '_temp_file_'


def copy_binary_file_to_hosts(ip_addresses):

  for ip in ip_addresses:
    command = ''
    command += ' scp -i ' + config.PRIVATE_KEY
    command += ' -o UserKnownHostsFile=/dev/null '
    command += ' -o StrictHostKeyChecking=no '
    command +=  config.SOURCE_PATH + 'Water '
    command += ' ubuntu@' + ip + ':' + config.REMOTE_PATH

    subprocess.call(command, shell=True)



def make_nodes_file_content(ip_addresses):

  string = ""
  for ip in ip_addresses:
    print ip
    string = string + ip + " cpu=8\n"

  file = open(temp_file_name, 'w+')
  file.write(string)
  file.close()


def copy_nodes_file_to_hosts(ip_addresses):
  make_nodes_file_content(ip_addresses)

  for ip in ip_addresses:
    command = ''
    command += ' scp -i ' + config.PRIVATE_KEY
    command += ' -o UserKnownHostsFile=/dev/null '
    command += ' -o StrictHostKeyChecking=no '
    command += temp_file_name
    command += ' ubuntu@' + ip + ':' + config.REMOTE_PATH + config.NODES_FILE_NAME

    subprocess.call(command, shell=True)

  subprocess.call(['rm', temp_file_name])


def run_experiment(ip):
  command = ''
  command += ' ssh -i ' + config.PRIVATE_KEY
  command += ' -o UserKnownHostsFile=/dev/null '
  command += ' -o StrictHostKeyChecking=no '
  command += ' ubuntu@' + ip
  command += ' \"cd ' + config.REMOTE_PATH + '; '
  command += ' mpirun -hostfile ' + config.NODES_FILE_NAME
  command += ' -np ' + str(config.INSTANCE_NUM)
  command += ' ./Water -scale ' + str(config.SCALE)
  command += ' -e ' + str(config.FRAME_NUM) + '\" '

  print command

  subprocess.call(command, shell=True)


def collect_output_data(ip_addresses):

  subprocess.call(['rm', '-rf', config.OUTPUT_NAME])
  subprocess.call(['mkdir', '-p', config.OUTPUT_NAME])

  process_num = 0
  for ip in ip_addresses:
    process_num += 1
    command = ''
    command += ' scp -r -i ' + config.PRIVATE_KEY
    command += ' -o UserKnownHostsFile=/dev/null '
    command += ' -o StrictHostKeyChecking=no '
    command += ' ubuntu@' + ip + ':' + config.REMOTE_PATH + config.OUTPUT_NAME + str(process_num)
    command += ' ' + config.OUTPUT_NAME

    subprocess.call(command,  shell=True)



