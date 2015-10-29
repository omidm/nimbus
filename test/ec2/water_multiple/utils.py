#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import sys
import os
import subprocess
import time

import config

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ec2

def run_experiment(controller_ip, controller_p_ip, worker_ips, worker_p_ips):

  assert(len(worker_ips) == config.WORKER_INSTANCE_NUM);
  assert(len(worker_p_ips) == config.WORKER_INSTANCE_NUM);

  worker_num = config.WORKER_INSTANCE_NUM * config.WORKER_PER_INSTANCE 
  run_controller(controller_ip, worker_num);
  time.sleep(5)

  for idx in range(0, config.WORKER_INSTANCE_NUM):
    time.sleep(config.WORKER_PER_INSTANCE - 1)
    ip = worker_ips[idx]
    p_ip = worker_p_ips[idx]
    for w in range(0, config.WORKER_PER_INSTANCE):
      run_worker(controller_p_ip, ip, p_ip, idx * config.WORKER_PER_INSTANCE + w);


def run_controller(controller_ip, worker_num):
  controller_command  = 'cd ' + config.EC2_NIMBUS_ROOT + config.REL_CONTROLLER_PATH + ';'
  controller_command += 'export DBG=' + config.DBG_MODE + ';'
  controller_command += 'export TTIMER=' + config.TTIMER_LEVEL + ';'
  controller_command += 'sudo ' + config.EC2_NIMBUS_ROOT + 'scripts/configure_tcp.sh;'
  controller_command += 'sudo sysctl -p;'
  controller_command += 'ulimit -c unlimited;'
  controller_command += './scheduler'
  controller_command += ' -p ' + str(config.FIRST_PORT)
  controller_command += ' -t ' + str(config.ASSIGNER_THREAD_NUM)
  controller_command += ' -a ' + str(config.BATCH_ASSIGN_NUM)
  controller_command += ' -c ' + str(config.COMMAND_BATCH_SIZE)
  controller_command += ' -w ' + str(worker_num)
  controller_command += ' --lb_period ' + str(config.LB_PERIOD)
  controller_command += ' --ft_period ' + str(config.FT_PERIOD)
  if config.ACTIVATE_LB:
    controller_command += ' --alb '
  if config.ACTIVATE_FT:
    controller_command += ' --aft '
  if config.DEACTIVATE_CONTROLLER_TEMPLATE:
    controller_command += ' --dct '
  if config.DEACTIVATE_COMPLEX_MEMOIZATION:
    controller_command += ' --dcm '
  if config.DEACTIVATE_BINDING_MEMOIZATION:
    controller_command += ' --dbm '
  if config.DEACTIVATE_WORKER_TEMPLATE:
    controller_command += ' --dwt '
  if config.DEACTIVATE_MEGA_RCR_JOB:
    controller_command += ' --dmr '
  if config.DEACTIVATE_DM_QUERY_CACHE:
    controller_command += ' --dqc '
  controller_command += ' &> ' + config.STD_OUT_LOG

  subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip, controller_command])

  print '** Controller Launched: ' + controller_ip


def run_worker(controller_p_ip, worker_ip, worker_p_ip, num):
  worker_command =  'cd ' + config.EC2_NIMBUS_ROOT + config.REL_WORKER_PATH + ';'
  worker_command += 'export DBG=' + config.DBG_MODE + ';'
  worker_command += 'export TTIMER=' + config.TTIMER_LEVEL + ';'
  worker_command += 'sudo ' + config.EC2_NIMBUS_ROOT + 'scripts/configure_tcp.sh;'
  worker_command += 'sudo sysctl -p;'
  worker_command += 'ulimit -c unlimited;'
  if config.RUN_WITH_TASKSET:
    worker_command += 'taskset -c ' + config.WORKER_TASKSET + ' ./worker'
  else:
    worker_command += './worker'
  worker_command += ' -p ' + str(config.FIRST_PORT + num + 1)
  worker_command += ' --ip ' + worker_p_ip
  worker_command += ' --cip ' + controller_p_ip
  worker_command += ' --cport ' + str(config.FIRST_PORT)
  worker_command += ' -s ' + str(config.SIMULATION_SCALE)
  worker_command += ' -e ' + str(config.FRAME_NUMBER)
  worker_command += ' --pnx ' + str(config.PART_X)
  worker_command += ' --pny ' + str(config.PART_Y)
  worker_command += ' --pnz ' + str(config.PART_Z)
  worker_command += ' --ppnx ' + str(config.PROJ_PART_X)
  worker_command += ' --ppny ' + str(config.PROJ_PART_Y)
  worker_command += ' --ppnz ' + str(config.PROJ_PART_Z)
  worker_command += ' --maxi ' + str(config.MAX_ITERATION)
  worker_command += ' --psl ' + str(config.PROJECTION_SMART_LEVEL)
  worker_command += ' --wl ' + str(config.WATER_LEVEL)
  worker_command += ' --ibatch ' + str(config.ITERATION_BATCH)
  worker_command += ' --othread ' + str(config.OTHREAD_NUM)
  if config.NO_EXECUTION_TEMPLATE:
    worker_command += ' --det '
  if not config.GLOBAL_WRITE:
    worker_command += ' --dgw '
  if config.NO_PROJ_BOTTLENECK:
    worker_command += ' --dpb '
  worker_command += ' &> ' + str(num + 1) + '_' + config.STD_OUT_LOG

  subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + worker_ip, worker_command])

  print '** Worker ' + str(num + 1) + ' Launched: ' + worker_ip


def terminate_experiment(controller_ip, worker_ips):

  terminate_controller(controller_ip);

  assert(len(worker_ips) == config.WORKER_INSTANCE_NUM);

  for ip in worker_ips:
    terminate_worker(ip);


def terminate_controller(controller_ip):
  controller_command =  'killall -v scheduler;'

  subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip, controller_command])

  print '** Controller Terminated: ' + controller_ip


def terminate_worker(worker_ip):
  worker_command =  'killall -v worker;'

  subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + worker_ip, worker_command])

  print '** Worker(s) Terminated: ' + worker_ip


def  test_nodes(node_ips):
  command  = 'cd ' + config.EC2_NIMBUS_ROOT + config.REL_WORKER_PATH + ';'
  command += 'pwd;'

  for ip in node_ips:
    subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip, command])

    print '** Node Tested: ' + ip


def collect_output_data(controller_ip, worker_ips):

  subprocess.call(['rm', '-rf', config.OUTPUT_PATH])
  subprocess.call(['mkdir', '-p', config.OUTPUT_PATH])

  subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip + ':' + config.EC2_NIMBUS_ROOT +
      config.REL_CONTROLLER_PATH + config.STD_OUT_LOG,
      config.OUTPUT_PATH])

  subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip + ':' + config.EC2_NIMBUS_ROOT +
      config.REL_CONTROLLER_PATH + config.LOAD_BALANCER_LOG,
      config.OUTPUT_PATH])

  subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip + ':' + config.EC2_NIMBUS_ROOT +
      config.REL_CONTROLLER_PATH + config.SCHED_PER_ITER_STAT_LOG,
      config.OUTPUT_PATH])

  for ip in worker_ips:
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


def clean_output_data(controller_ip, worker_ips):
  controller_path = config.EC2_NIMBUS_ROOT + config.REL_CONTROLLER_PATH;
  controller_command  =  ''
  controller_command +=  'rm -rf ' + controller_path + config.STD_OUT_LOG + ';'
  controller_command +=  'rm -rf ' + controller_path + '*.txt;'
  controller_command +=  'rm -rf ' + controller_path + '*log*;'
  controller_command +=  'rm -rf ' + controller_path + 'core;'

  subprocess.Popen(['ssh', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip, controller_command])

  print '** Controller Cleaned: ' + controller_ip

  for ip in worker_ips:
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
  
    print '** Worker Cleaned: ' + ip









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
  command += 'cd ' + config.REL_CONTROLLER_PATH + ';'
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
        'ubuntu@' + source_ip + ':' + config.EC2_NIMBUS_ROOT + config.REL_CONTROLLER_PATH + config.CONTROLLER_BINARY,
        'ubuntu@' + dest_ip + ':' + config.EC2_NIMBUS_ROOT + config.REL_CONTROLLER_PATH])

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
#       temp_folder_name + config.REL_CONTROLLER_PATH])
#   subprocess.call(['mkdir', '-p',
#       temp_folder_name + config.REL_WORKER_PATH])
#   subprocess.call(['cp',
#       config.SOURCE_NIMBUS_ROOT + config.NIMBUS_LIB,
#       temp_folder_name])
#   subprocess.call(['cp',
#       config.SOURCE_NIMBUS_ROOT  + config.REL_APPLICATION_PATH + config.APPLICATION_LIB,
#       temp_folder_name + config.REL_APPLICATION_PATH])
#   subprocess.call(['cp',
#       config.SOURCE_NIMBUS_ROOT  + config.REL_CONTROLLER_PATH + config.CONTROLLER_BINARY,
#       temp_folder_name + config.REL_CONTROLLER_PATH])
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






