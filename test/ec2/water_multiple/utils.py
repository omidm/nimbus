#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import sys
import os
import subprocess
import time

import config

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ec2


def build_binaries(ip_address):

# Thing you need to do first on all hosts:
#     git pull
#     Build PhysBAM library
#
# Thing you need to do on source host:
#  In application/water_multiple:
#     generate the data_def and region_def files
#     update parameter.h
#     create the make file with ccmake

  command = ''
  command += 'cd ' + config.EC2_NIMBUS_ROOT + ';'
  command += 'make clean-hard;'
  command += 'make;'
  command += 'cd ' + config.REL_APPLICATION_PATH + ';'
  command += 'make clean;'
  command += 'make -j 12;'
  command += 'cd -;'
  command += 'cd ' + config.REL_WORKER_PATH + ';'
  command += 'make clean;'
  command += 'make -j 12;'
  command += 'cd -;'
  command += 'cd ' + config.REL_SCHEDULER_PATH + ';'
  command += 'make clean;'
  command += 'make -j 12;'
  command += 'cd -;'

  subprocess.call(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + ip_address, command])


def distribute_binaries(source_ip, dest_ips):
  for dest_ip in dest_ips:
    subprocess.call(['scp', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + source_ip + ':' + config.EC2_NIMBUS_ROOT + config.NIMBUS_LIB,
        'ubuntu@' + dest_ip + ':' + config.EC2_NIMBUS_ROOT])

  for dest_ip in dest_ips:
    subprocess.call(['scp', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + source_ip + ':' + config.EC2_NIMBUS_ROOT + config.REL_APPLICATION_PATH + config.APPLICATION_LIB,
        'ubuntu@' + dest_ip + ':' + config.EC2_NIMBUS_ROOT + config.REL_APPLICATION_PATH])

  for dest_ip in dest_ips:
    subprocess.call(['scp', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + source_ip + ':' + config.EC2_NIMBUS_ROOT + config.REL_SCHEDULER_PATH + config.SCHEDULER_BINARY,
        'ubuntu@' + dest_ip + ':' + config.EC2_NIMBUS_ROOT + config.REL_SCHEDULER_PATH])

  for dest_ip in dest_ips:
    subprocess.call(['scp', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + source_ip + ':' + config.EC2_NIMBUS_ROOT + config.REL_WORKER_PATH + config.WORKER_BINARY,
        'ubuntu@' + dest_ip + ':' + config.EC2_NIMBUS_ROOT + config.REL_WORKER_PATH])

def run_experiment(scheduler_ip, worker_ips):

  worker_num = len(worker_ips)

  run_scheduler(scheduler_ip, worker_num);
  time.sleep(10)

  num = 0;
  for ip in worker_ips:
    num += 1
    run_worker(scheduler_ip, ip, num);


def run_scheduler(scheduler_ip, worker_num):

  scheduler_command =  'cd ' + config.EC2_NIMBUS_ROOT + config.REL_SCHEDULER_PATH + ';'
  scheduler_command += 'export DBG=error;'
  scheduler_command += 'sudo ' + config.EC2_NIMBUS_ROOT + 'scripts/configure_tcp.sh;'
  scheduler_command += 'ulimit -c unlimited;'
  scheduler_command += './scheduler'
  scheduler_command += ' -port ' + str(config.FIRST_PORT)
  scheduler_command += ' -wn ' + str(worker_num)
  scheduler_command += ' -tn ' + str(config.ASSIGNER_THREAD_NUM)
  scheduler_command += ' -an ' + str(config.BATCH_ASSIGN_NUM)
  scheduler_command += ' &> ' + config.LOG_FILE_NAME

  subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip, scheduler_command])

  print '** Scheduler Launched: ' + scheduler_ip


def run_worker(scheduler_ip, worker_ip, num):
  worker_command =  'cd ' + config.EC2_NIMBUS_ROOT + config.REL_WORKER_PATH + ';'
  worker_command += 'export DBG=error;'
  worker_command += 'sudo ' + config.EC2_NIMBUS_ROOT + 'scripts/configure_tcp.sh;'
  worker_command += 'ulimit -c unlimited;'
  worker_command += './worker'
  worker_command += ' -port ' + str(config.FIRST_PORT + num)
  worker_command += ' -ip ' + worker_ip
  worker_command += ' -sip ' + scheduler_ip
  worker_command += ' -sport ' + str(config.FIRST_PORT)
  worker_command += ' -port ' + str(config.FIRST_PORT + num)
  worker_command += ' -othread ' + str(config.OTHREAD_NUM)
  worker_command += ' &> ' + str(num) + '_' + config.LOG_FILE_NAME

  subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + worker_ip, worker_command])

  print '** Worker ' + str(num) + ' Launched: ' + worker_ip

def terminate_experiment(scheduler_ip, worker_ips):

  worker_num = len(worker_ips)

  terminate_scheduler(scheduler_ip);

  num = 0;
  for ip in worker_ips:
    num += 1
    terminate_worker(ip, num);


def terminate_scheduler(scheduler_ip):
  scheduler_command =  'killall -v scheduler;'

  subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip, scheduler_command])

  print '** Scheduler Terminated: ' + scheduler_ip


def terminate_worker(worker_ip, num):
  worker_command =  'killall -v worker;'

  subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + worker_ip, worker_command])

  print '** Worker ' + str(num) + ' Terminated: ' + worker_ip

def  test_workers(worker_ips):
  worker_command =  'cd ' + config.EC2_NIMBUS_ROOT + config.REL_WORKER_PATH + ';'
  worker_command += 'pwd;'

  num = 0;
  for ip in worker_ips:
    num += 1;
    subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip, worker_command])

    print '** Worker ' + str(num) + ' Tested: ' + ip





def collect_output_data(scheduler_ip, worker_ips):

  subprocess.call(['rm', '-rf', config.OUTPUT_PATH])
  subprocess.call(['mkdir', '-p', config.OUTPUT_PATH])

  subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip + ':' + config.EC2_NIMBUS_ROOT +
      config.REL_SCHEDULER_PATH + config.LOG_FILE_NAME,
      config.OUTPUT_PATH])

  subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip + ':' + config.EC2_NIMBUS_ROOT +
      config.REL_SCHEDULER_PATH + config.SCHED_LOG_NAME_1,
      config.OUTPUT_PATH])

  subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip + ':' + config.EC2_NIMBUS_ROOT +
      config.REL_SCHEDULER_PATH + config.SCHED_LOG_NAME_2,
      config.OUTPUT_PATH])

  subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip + ':' + config.EC2_NIMBUS_ROOT +
      config.REL_SCHEDULER_PATH + config.SCHED_LOG_NAME_3,
      config.OUTPUT_PATH])

  subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip + ':' + config.EC2_NIMBUS_ROOT +
      config.REL_SCHEDULER_PATH + 'core',
      config.OUTPUT_PATH])

  num = 0
  for ip in worker_ips:
    num += 1
    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        # config.REL_WORKER_PATH + str(num) + '_' + config.LOG_FILE_NAME,
        config.REL_WORKER_PATH  + '*_' + config.LOG_FILE_NAME,
        config.OUTPUT_PATH])

    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        # config.REL_WORKER_PATH + config.WORKER_LOG_FILE_NAME + str(config.FIRST_PORT + num),
        config.REL_WORKER_PATH + config.WORKER_LOG_FILE_NAME + '*',
        config.OUTPUT_PATH])

    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        config.REL_WORKER_PATH + 'event_fe.txt',
        config.OUTPUT_PATH + str(num) + '_event_fe.txt'])

    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        config.REL_WORKER_PATH + 'event_be.txt',
        config.OUTPUT_PATH + str(num) + '_event_be.txt'])

    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        config.REL_WORKER_PATH + 'core',
        config.OUTPUT_PATH + str(num) + '_core'])

    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        config.REL_WORKER_PATH + 'cache_objects.txt',
        config.OUTPUT_PATH + str(num) + '_cache_objects.txt'])

    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        config.REL_WORKER_PATH + 'data_objects.txt',
        config.OUTPUT_PATH + str(num) + '_data_objects.txt'])

    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        config.REL_WORKER_PATH + 'cache_behavior.txt',
        config.OUTPUT_PATH + str(num) + '_cache_behavior.txt'])



def clean_output_data(scheduler_ip, worker_ips):
  scheduler_path = config.EC2_NIMBUS_ROOT + config.REL_SCHEDULER_PATH;
  scheduler_command  =  ''
  scheduler_command +=  'rm -rf ' + scheduler_path + config.LOG_FILE_NAME + ';'
  scheduler_command +=  'rm -rf ' + scheduler_path + '*log*;'
  scheduler_command +=  'rm -rf ' + scheduler_path + 'core;'

  subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip, scheduler_command])

  print '** Scheduler Cleaned: ' + scheduler_ip

  num = 0
  for ip in worker_ips:
    num += 1
    worker_path = config.EC2_NIMBUS_ROOT + config.REL_WORKER_PATH;
    worker_command  =  ''
    worker_command +=  'rm -rf ' + worker_path + '*_' + config.LOG_FILE_NAME + ';'
    worker_command +=  'rm -rf ' + worker_path + config.WORKER_LOG_FILE_NAME + '*;'
    worker_command +=  'rm -rf ' + worker_path + 'event_*;'
    worker_command +=  'rm -rf ' + worker_path + 'core;'
  
    subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip, worker_command])
  
    print '** Worker ' + str(num) + ' Cleaned: ' + ip

# ssh -i omidm-sing-key-pair-us-west-2.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ubuntu@<ip>

# def create_folder_of_binaries():
#   temp_folder_name = '_temp_folder_/'
# 
#   subprocess.call(['rm', '-rf',
#       temp_folder_name])
#   subprocess.call(['mkdir', '-p',
#       temp_folder_name])
#   subprocess.call(['mkdir', '-p',
#       temp_folder_name + config.REL_PHYSBAM_PATH])
#   subprocess.call(['mkdir', '-p',
#       temp_folder_name + config.REL_PHYSBAM_PATH2])
#   subprocess.call(['mkdir', '-p',
#       temp_folder_name + config.REL_APPLICATION_PATH])
#   subprocess.call(['mkdir', '-p',
#       temp_folder_name + config.REL_SCHEDULER_PATH])
#   subprocess.call(['mkdir', '-p',
#       temp_folder_name + config.REL_WORKER_PATH])
#   subprocess.call(['cp',
#       config.SOURCE_NIMBUS_ROOT + config.NIMBUS_LIB,
#       temp_folder_name])
#   subprocess.call(['cp',
#       config.SOURCE_NIMBUS_ROOT  + config.REL_APPLICATION_PATH + config.APPLICATION_LIB,
#       temp_folder_name + config.REL_APPLICATION_PATH])
#   subprocess.call(['cp',
#       config.SOURCE_NIMBUS_ROOT  + config.REL_SCHEDULER_PATH + config.SCHEDULER_BINARY,
#       temp_folder_name + config.REL_SCHEDULER_PATH])
#   subprocess.call(['cp',
#       config.SOURCE_NIMBUS_ROOT  + config.REL_WORKER_PATH + config.WORKER_BINARY,
#       temp_folder_name + config.REL_WORKER_PATH])
# 
#   string  = 'cp '
#   string += config.SOURCE_NIMBUS_ROOT  + config.REL_PHYSBAM_PATH + config.PHYSBAM_LIB
#   string += ' ' + temp_folder_name + config.REL_PHYSBAM_PATH
#   subprocess.call(string, shell=True)
# 
#   string  = 'cp '
#   string += config.SOURCE_NIMBUS_ROOT  + config.REL_PHYSBAM_PATH2 + config.PHYSBAM_LIB
#   string += ' ' + temp_folder_name + config.REL_PHYSBAM_PATH2
#   subprocess.call(string, shell=True)
# 
#   return temp_folder_name
# 
# 
# def copy_binary_folder_to_hosts(ip_addresses):
#   temp_folder = create_folder_of_binaries()
# 
#   for ip in ip_addresses:
#     subprocess.call(['scp', '-r', '-i', config.PRIVATE_KEY,
#         '-o', 'UserKnownHostsFile=/dev/null',
#         '-o', 'StrictHostKeyChecking=no',
#         temp_folder,
#         'ubuntu@' + ip + ':~/' + config.EC2_NIMBUS_ROOT])
# 
#   subprocess.call(['rm', '-rf', temp_folder])






