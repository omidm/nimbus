#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import sys
import os
import subprocess

import config

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ec2



def create_folder_of_binaries():
  temp_folder_name = '_temp_folder_/'

  subprocess.call(['rm', '-rf',
      temp_folder_name])
  subprocess.call(['mkdir', '-p',
      temp_folder_name])
  subprocess.call(['mkdir', '-p',
      temp_folder_name + config.REL_APPLICATION_PATH])
  subprocess.call(['mkdir', '-p',
      temp_folder_name + config.REL_SCHEDULER_PATH])
  subprocess.call(['mkdir', '-p',
      temp_folder_name + config.REL_WORKER_PATH])
  subprocess.call(['cp',
      config.SOURCE_NIMBUS_ROOT + config.NIMBUS_LIB,
      temp_folder_name])
  subprocess.call(['cp',
      config.SOURCE_NIMBUS_ROOT  + config.REL_APPLICATION_PATH + config.APPLICATION_LIB,
      temp_folder_name + config.REL_APPLICATION_PATH])
  subprocess.call(['cp',
      config.SOURCE_NIMBUS_ROOT  + config.REL_SCHEDULER_PATH + config.SCHEDULER_BINARY,
      temp_folder_name + config.REL_SCHEDULER_PATH])
  subprocess.call(['cp',
      config.SOURCE_NIMBUS_ROOT  + config.REL_WORKER_PATH + config.WORKER_BINARY,
      temp_folder_name + config.REL_WORKER_PATH])

  return temp_folder_name


def copy_binary_folder_to_hosts(ip_addresses):

  temp_folder = create_folder_of_binaries()

  for ip in ip_addresses:
    subprocess.call(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        temp_folder,
        'ubuntu@' + ip + ':~/' + config.EC2_FOLDER_NAME])

  subprocess.call(['rm', '-rf', temp_folder])



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
    subprocess.call(['scp', '-i', '/home/omidm/.ssh/' + config.KEY_NAME + '.pem',
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        temp_file_name, 'ubuntu@' + ip + ':~/' + physbam_config.DIRECTORY_PATH +
        physbam_config.NODES_FILE_NAME])

  subprocess.call(['rm', temp_file_name])


def run_experiment(ip):
  print "omid"


  command =  'cd ' + physbam_config.DIRECTORY_PATH + ';'
  command += 'mpirun -hostfile ' + physbam_config.NODES_FILE_NAME
  command += ' -np ' + str(config.INSTANCE_NUM)
  command += ' ./Water -scale 128 -e 40'

  subprocess.call(['ssh', '-i', '/home/omidm/.ssh/' + config.KEY_NAME + '.pem',
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + ip, command])


def collect_output_data(ip_addresses):

  subprocess.call(['rm', '-rf', physbam_config.OUTPUT_PATH])
  subprocess.call(['mkdir', '-p', physbam_config.OUTPUT_PATH])

  process_num = 0
  for ip in ip_addresses:
    process_num += 1
    subprocess.call(['scp', '-r', '-i', '/home/omidm/.ssh/' + config.KEY_NAME + '.pem',
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':~/' + physbam_config.DIRECTORY_PATH +
        physbam_config.OUTPUT_PATH + str(process_num),
        physbam_config.OUTPUT_PATH])



