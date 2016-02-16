#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import sys
import os
import subprocess
import time

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



def run_experiment(scheduler_ip, worker_ips):

  worker_num = len(worker_ips)

  print '** Scheduler'
  print scheduler_ip
  print '***'
  scheduler_command =  'cd ' + config.EC2_FOLDER_NAME + config.REL_SCHEDULER_PATH + ';'
  scheduler_command += './scheduler'
  scheduler_command += ' -port ' + str(config.FIRST_PORT)
  scheduler_command += ' -wn ' + str(worker_num)
  scheduler_command += ' > ' + config.LOG_FILE_NAME

  subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip, scheduler_command])

  time.sleep(10)

  num = 0;
  for ip in worker_ips:
    time.sleep(5)
    print '** Worker'
    print ip
    print '***'
    num += 1
    worker_command =  'cd ' + config.EC2_FOLDER_NAME + config.REL_WORKER_PATH + ';'
    worker_command += './worker'
    worker_command += ' -port ' + str(config.FIRST_PORT + num)
    worker_command += ' -ip ' + ip
    worker_command += ' -sip ' + scheduler_ip
    worker_command += ' -sport ' + str(config.FIRST_PORT)
    worker_command += './worker -port ' + str(config.FIRST_PORT + num)
    worker_command += ' > ' + str(num) + '_' + config.LOG_FILE_NAME

 
    if num == len(worker_ips):
      subprocess.call(['ssh', '-i', config.PRIVATE_KEY,
          '-o', 'UserKnownHostsFile=/dev/null',
          '-o', 'StrictHostKeyChecking=no',
          'ubuntu@' + ip, worker_command])
    else:
      subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
          '-o', 'UserKnownHostsFile=/dev/null',
          '-o', 'StrictHostKeyChecking=no',
          'ubuntu@' + ip, worker_command])


def collect_output_data(scheduler_ip, worker_ips):

  subprocess.call(['rm', '-rf', config.OUTPUT_PATH])
  subprocess.call(['mkdir', '-p', config.OUTPUT_PATH])

  subprocess.call(['scp', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip + ':~/' + config.EC2_FOLDER_NAME +
      config.REL_SCHEDULER_PATH + config.LOG_FILE_NAME,
      config.OUTPUT_PATH])

  num = 0
  for ip in worker_ips:
    num += 1
    subprocess.call(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':~/' + config.EC2_FOLDER_NAME +
        config.REL_WORKER_PATH + str(num) + '_' + config.LOG_FILE_NAME,
        config.OUTPUT_PATH])



