#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import sys
import os
import subprocess
import time

import config

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ec2

def run_experiment(scheduler_ip, scheduler_p_ip, worker_ips, worker_p_ips):

  worker_num = len(worker_ips)

  run_scheduler(scheduler_ip, worker_num);
  time.sleep(10)

  idx = 0;
  for idx in range(0, len(worker_ips)):
    ip = worker_ips[idx]
    p_ip = worker_p_ips[idx]
    run_worker(scheduler_p_ip, ip, p_ip, idx+1);


def run_scheduler(scheduler_ip, worker_num):
  assert(worker_num == config.WORKER_NUM)

  scheduler_command =  'cd ' + config.EC2_NIMBUS_ROOT + config.REL_SCHEDULER_PATH + ';'
  scheduler_command += 'export DBG=error;'
  scheduler_command += 'sudo ' + config.EC2_NIMBUS_ROOT + 'scripts/configure_tcp.sh;'
  scheduler_command += 'sudo sysctl -p;'
  scheduler_command += 'ulimit -c unlimited;'
  scheduler_command += './scheduler'
  scheduler_command += ' -p ' + str(config.FIRST_PORT)
  scheduler_command += ' -w ' + str(config.WORKER_NUM)
  scheduler_command += ' -t ' + str(config.ASSIGNER_THREAD_NUM)
  scheduler_command += ' -a ' + str(config.BATCH_ASSIGN_NUM)
  scheduler_command += ' -c ' + str(config.COMMAND_BATCH_SIZE)
  if config.ACTIVATE_LB:
    scheduler_command += ' --alb '
  if config.ACTIVATE_FT:
    scheduler_command += ' --aft '
  scheduler_command += ' --lb_period ' + str(config.LB_PERIOD)
  scheduler_command += ' --ft_period ' + str(config.FT_PERIOD)
  if config.DEACTIVATE_CONTROLLER_TEMPLATE:
    scheduler_command += ' --dct '
  if config.DEACTIVATE_COMPLEX_MEMOIZATION:
    scheduler_command += ' --dcm '
  if config.DEACTIVATE_BINDING_MEMOIZATION:
    scheduler_command += ' --dbm '
  if config.DEACTIVATE_WORKER_TEMPLATE:
    scheduler_command += ' --dwt '
  if config.DEACTIVATE_MEGA_RCR_JOB:
    scheduler_command += ' --dmr '
  if config.DEACTIVATE_DM_QUERY_CACHE:
    scheduler_command += ' --dqc '
  scheduler_command += ' &> ' + config.STD_OUT_LOG

  subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip, scheduler_command])

  print '** Scheduler Launched: ' + scheduler_ip


def run_worker(scheduler_p_ip, worker_ip, worker_p_ip, num):
  worker_command =  'cd ' + config.EC2_NIMBUS_ROOT + config.REL_WORKER_PATH + ';'
  worker_command += 'export DBG=error;'
  worker_command += 'sudo ' + config.EC2_NIMBUS_ROOT + 'scripts/configure_tcp.sh;'
  worker_command += 'sudo sysctl -p;'
  worker_command += 'ulimit -c unlimited;'
  if config.RUN_WITH_TASKSET:
    worker_command += 'taskset -c ' + config.WORKER_TASKSET + ' ./worker'
  else:
    worker_command += './worker'
  worker_command += ' -s ' + str(config.SIMULATION_SCALE)
  worker_command += ' -e ' + str(config.FRAME_NUMBER)
  worker_command += ' --pnx ' + str(config.PART_X)
  worker_command += ' --pny ' + str(config.PART_Y)
  worker_command += ' --pnz ' + str(config.PART_Z)
  worker_command += ' --ppnx ' + str(config.PROJ_PART_X)
  worker_command += ' --ppny ' + str(config.PROJ_PART_Y)
  worker_command += ' --ppnz ' + str(config.PROJ_PART_Z)
  worker_command += ' --maxi ' + str(config.MAX_ITERATION)
  worker_command += ' --ibatch ' + str(config.ITERATION_BATCH)
  if config.WRITE_PER_PART:
    worker_command += ' --dgw '
  if config.NO_PROJ_BOTTLENECK:
    worker_command += ' --dpb '
  worker_command += ' -p ' + str(config.FIRST_PORT + num)
  worker_command += ' --ip ' + worker_p_ip
  worker_command += ' --cip ' + scheduler_p_ip
  worker_command += ' --cport ' + str(config.FIRST_PORT)
  worker_command += ' --othread ' + str(config.OTHREAD_NUM)
  worker_command += ' &> ' + str(num) + '_' + config.STD_OUT_LOG

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

def  test_nodes(node_ips):
  worker_command =  'cd ' + config.EC2_NIMBUS_ROOT + config.REL_WORKER_PATH + ';'
  worker_command += 'pwd;'

  num = 0;
  for ip in node_ips:
    num += 1;
    subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip, worker_command])

    print '** Node ' + str(num) + ' Tested: ' + ip





def collect_output_data(scheduler_ip, worker_ips):

  subprocess.call(['rm', '-rf', config.OUTPUT_PATH])
  subprocess.call(['mkdir', '-p', config.OUTPUT_PATH])

  subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip + ':' + config.EC2_NIMBUS_ROOT +
      config.REL_SCHEDULER_PATH + config.STD_OUT_LOG,
      config.OUTPUT_PATH])

  subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip + ':' + config.EC2_NIMBUS_ROOT +
      config.REL_SCHEDULER_PATH + config.LOAD_BALANCER_LOG,
      config.OUTPUT_PATH])

  subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + scheduler_ip + ':' + config.EC2_NIMBUS_ROOT +
      config.REL_SCHEDULER_PATH + config.SCHED_PER_ITER_STAT_LOG,
      config.OUTPUT_PATH])

  num = 0
  for ip in worker_ips:
    num += 1
    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        config.REL_WORKER_PATH  + '*_' + config.STD_OUT_LOG,
        config.OUTPUT_PATH])

    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        config.REL_WORKER_PATH + '*_' + config.WORKER_LB_LOG,
        config.OUTPUT_PATH])

    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        config.REL_WORKER_PATH + '*_' + config.WORKER_FINAL_STAT_LOG,
        config.OUTPUT_PATH])

    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        config.REL_WORKER_PATH + '*_' + config.WORKER_PER_ITER_STAT_LOG,
        config.OUTPUT_PATH])

    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        config.REL_WORKER_PATH + 'mpi*.log',
        config.OUTPUT_PATH])


#    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
#        '-o', 'UserKnownHostsFile=/dev/null',
#        '-o', 'StrictHostKeyChecking=no',
#        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
#        config.REL_WORKER_PATH + 'core',
#        config.OUTPUT_PATH + str(num) + '_core'])
#
#    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
#        '-o', 'UserKnownHostsFile=/dev/null',
#        '-o', 'StrictHostKeyChecking=no',
#        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
#        config.REL_WORKER_PATH + '*_cache_time.txt',
#        config.OUTPUT_PATH])
#
#    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
#        '-o', 'UserKnownHostsFile=/dev/null',
#        '-o', 'StrictHostKeyChecking=no',
#        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
#        config.REL_WORKER_PATH + '*_cache_objects.txt',
#        config.OUTPUT_PATH])
#
#    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
#        '-o', 'UserKnownHostsFile=/dev/null',
#        '-o', 'StrictHostKeyChecking=no',
#        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
#        config.REL_WORKER_PATH + '*_cache_behavior.txt',
#        config.OUTPUT_PATH])
#
#    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
#        '-o', 'UserKnownHostsFile=/dev/null',
#        '-o', 'StrictHostKeyChecking=no',
#        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
#        config.REL_WORKER_PATH  + 'log_wdx*',
#        config.OUTPUT_PATH])


def clean_output_data(scheduler_ip, worker_ips):
  scheduler_path = config.EC2_NIMBUS_ROOT + config.REL_SCHEDULER_PATH;
  scheduler_command  =  ''
  scheduler_command +=  'rm -rf ' + scheduler_path + config.STD_OUT_LOG + ';'
  scheduler_command +=  'rm -rf ' + scheduler_path + '*.txt;'
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
    worker_command +=  'rm -rf ' + worker_path + '*_' + config.STD_OUT_LOG + ';'
    worker_command +=  'rm -rf ' + worker_path + '*.txt;'
    worker_command +=  'rm -rf ' + worker_path + '*log*;'
    worker_command +=  'rm -rf ' + worker_path + 'core;'
    worker_command +=  'rm -rf ' + worker_path + 'split_output/;'
    worker_command +=  'rm -rf ' + worker_path + 'output/;'
    worker_command +=  'rm -rf ' + worker_path + '_db*;'
  
    subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip, worker_command])
  
    print '** Worker ' + str(num) + ' Cleaned: ' + ip


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






