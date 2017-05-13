#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import sys
import os
import subprocess
import time

import config
import ec2

# Logging configurations
STD_OUT_LOG                     = 'ec2_log.txt'
LOAD_BALANCER_LOG               = 'load_balancer_log'
CONTROLLER_PER_ITER_STAT_LOG    = 'controller_stats.txt'
WORKER_FINAL_STAT_LOG           = 'time_per_thread.txt' 
WORKER_PER_ITER_STAT_LOG        = 'main_timers.txt' 
OUTPUT_PATH                     = 'output/'

# Path configuration
NIMBUS_ROOT                     = '~/nimbus/'
REL_LOGS_PATH                   = 'logs/'
REL_CONTROLLER_PATH             = 'nodes/nimbus_controller/'
REL_WORKER_PATH                 = 'nodes/nimbus_worker/'


# Forming the worker arguments based on options

# lr
LR_REL_APP_PATH = 'applications/ml/logistic_regression/liblr.so' 
LR_APP_OPTIONS  = ' '
LR_APP_OPTIONS += ' -d ' + str(config.DIMENSION)
LR_APP_OPTIONS += ' -i ' + str(config.ITERATION_NUM)
LR_APP_OPTIONS += ' -s ' + str(config.SAMPLE_NUM_M)
LR_APP_OPTIONS += ' -p ' + str(config.PARTITION_NUM)
LR_APP_OPTIONS += ' -w ' + str(config.SPIN_WAIT_US)
LR_APP_OPTIONS += ' -r ' + str(config.REDUCTION_PARTITION_NUM)
if config.DEACTIVATE_AUTOMATIC_REDUCTION:
  LR_APP_OPTIONS += ' --dar '
if config.DEACTIVATE_REDUCTION_COMBINER:
  LR_APP_OPTIONS += ' --drc '

# k-means
KM_REL_APP_PATH = 'applications/ml/k_means/libk_means.so' 
KM_APP_OPTIONS  = ' '
KM_APP_OPTIONS += ' -d ' + str(config.DIMENSION)
KM_APP_OPTIONS += ' -c ' + str(config.CLUSTER_NUM)
KM_APP_OPTIONS += ' -i ' + str(config.ITERATION_NUM)
KM_APP_OPTIONS += ' -s ' + str(config.SAMPLE_NUM_M)
KM_APP_OPTIONS += ' -p ' + str(config.PARTITION_NUM)
KM_APP_OPTIONS += ' -w ' + str(config.SPIN_WAIT_US)
KM_APP_OPTIONS += ' -r ' + str(config.REDUCTION_PARTITION_NUM)
if config.DEACTIVATE_AUTOMATIC_REDUCTION:
  KM_APP_OPTIONS += ' --dar '
if config.DEACTIVATE_REDUCTION_COMBINER:
  KM_APP_OPTIONS += ' --drc '

# water
WATER_REL_APP_PATH = 'applications/physbam/water/libwater_app.so'
WATER_APP_OPTIONS  = ' '
WATER_APP_OPTIONS += ' -s ' + str(config.SIMULATION_SCALE)
WATER_APP_OPTIONS += ' -e ' + str(config.FRAME_NUMBER)
WATER_APP_OPTIONS += ' --pnx ' + str(config.PART_X)
WATER_APP_OPTIONS += ' --pny ' + str(config.PART_Y)
WATER_APP_OPTIONS += ' --pnz ' + str(config.PART_Z)
WATER_APP_OPTIONS += ' --ppnx ' + str(config.PROJ_PART_X)
WATER_APP_OPTIONS += ' --ppny ' + str(config.PROJ_PART_Y)
WATER_APP_OPTIONS += ' --ppnz ' + str(config.PROJ_PART_Z)
WATER_APP_OPTIONS += ' --maxi ' + str(config.MAX_ITERATION)
WATER_APP_OPTIONS += ' --psl ' + str(config.PROJECTION_SMART_LEVEL)
WATER_APP_OPTIONS += ' --wl ' + str(config.WATER_LEVEL)
WATER_APP_OPTIONS += ' --ibatch ' + str(config.ITERATION_BATCH)
if not config.GLOBAL_WRITE:
  WATER_APP_OPTIONS += ' --dgw '
if config.NO_PROJ_BOTTLENECK:
  WATER_APP_OPTIONS += ' --dpb '

# heat
HEAT_REL_APP_PATH = 'applications/stencilprobe/heat/libheat.so'
HEAT_APP_OPTIONS  = ' '
HEAT_APP_OPTIONS += ' -i ' + str(config.ITER_NUM)
HEAT_APP_OPTIONS += ' -x ' + str(config.NX)
HEAT_APP_OPTIONS += ' -y ' + str(config.NY)
HEAT_APP_OPTIONS += ' -z ' + str(config.NZ)
HEAT_APP_OPTIONS += ' --pnx ' + str(config.PNX)
HEAT_APP_OPTIONS += ' --pny ' + str(config.PNY)
HEAT_APP_OPTIONS += ' --pnz ' + str(config.PNZ)
HEAT_APP_OPTIONS += ' -w ' + str(config.SW_US)
HEAT_APP_OPTIONS += ' -b ' + str(config.BW)


if (config.APPLICATION == 'lr'):
  REL_APP_PATH = LR_REL_APP_PATH
  APP_OPTIONS  = LR_APP_OPTIONS
elif (config.APPLICATION == 'k-means'):
  REL_APP_PATH = KM_REL_APP_PATH
  APP_OPTIONS  = KM_APP_OPTIONS
elif (config.APPLICATION == 'water'):
  REL_APP_PATH = WATER_REL_APP_PATH
  APP_OPTIONS  = WATER_APP_OPTIONS
elif (config.APPLICATION == 'heat'):
  REL_APP_PATH = HEAT_REL_APP_PATH
  APP_OPTIONS  = HEAT_APP_OPTIONS
else:
  print "ERROR: unknown application: " + config.APPLICATION
  exit(1)


def start_experiment(controller_ip, controller_p_ip, worker_ips, worker_p_ips):
  start_controller(controller_ip, config.WORKER_NUM);

  assert(config.WORKER_NUM <= len(worker_ips))
  for idx in range(0, config.WORKER_NUM):
    ip = worker_ips[idx]
    p_ip = worker_p_ips[idx]
    start_worker(controller_p_ip, ip, p_ip, idx+1);


def start_controller(controller_ip, worker_num):
  assert(worker_num == config.WORKER_NUM)

  controller_command =  'cd ' + NIMBUS_ROOT + ';'
  controller_command += 'export DBG=' + config.DBG_MODE + ';'
  controller_command += 'export TTIMER=' + config.TTIMER_LEVEL + ';'
  controller_command += 'sudo scripts/configure_tcp.sh &> ' + STD_OUT_LOG + ';'
  controller_command += 'sudo sysctl -p &>> ' + STD_OUT_LOG + ';'
  controller_command += 'ulimit -c unlimited;'
  controller_command += 'scripts/start-controller.sh'
  controller_command += ' -p ' + str(config.FIRST_PORT)
  controller_command += ' -w ' + str(config.WORKER_NUM)
  controller_command += ' -t ' + str(config.ASSIGNER_THREAD_NUM)
  controller_command += ' -a ' + str(config.BATCH_ASSIGN_NUM)
  controller_command += ' -c ' + str(config.COMMAND_BATCH_SIZE)
  controller_command += ' --split ' + config.SPLIT_ARGS
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
  if config.DEACTIVATE_CASCADED_BINDING:
    controller_command += ' --dcb '
  controller_command += ' &>> ' + STD_OUT_LOG

  # print '** Starting controller: ' + controller_ip
  subprocess.Popen(['ssh', '-q', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip, controller_command])


def start_worker(controller_p_ip, worker_ip, worker_p_ip, num):
  worker_command =  'cd ' + NIMBUS_ROOT + ';'
  worker_command += 'export DBG=' + config.DBG_MODE + ';'
  worker_command += 'export TTIMER=' + config.TTIMER_LEVEL + ';'
  worker_command += 'sudo scripts/configure_tcp.sh &> ' + str(num) + '_' + STD_OUT_LOG + ';'
  worker_command += 'sudo sysctl -p &>>' + str(num) + '_' + STD_OUT_LOG + ';'
  worker_command += 'ulimit -c unlimited;'
  if config.DISABLE_HYPERTHREADING:
    worker_command += 'sudo scripts/disable_hyperthreading.sh &>> ' + str(num) + '_' + STD_OUT_LOG + ';'
  if config.RUN_WITH_TASKSET:
    worker_command += 'taskset -c ' + config.WORKER_TASKSET + ' ./scripts/start-workers.sh 1'
  else:
    worker_command += './scripts/start-workers.sh 1'

  worker_command += ' -p ' + str(config.FIRST_PORT + num)
  worker_command += ' --ip ' + worker_p_ip
  worker_command += ' --cip ' + controller_p_ip
  worker_command += ' --cport ' + str(config.FIRST_PORT)
  worker_command += ' --othread ' + str(config.OTHREAD_NUM)
  if config.DEACTIVATE_EXECUTION_TEMPLATE:
    worker_command += ' --det '
  worker_command += ' -l ' + REL_APP_PATH
  worker_command += ' ' + APP_OPTIONS
  worker_command += ' &>> ' + str(num) + '_' + STD_OUT_LOG

  # print '** Starting worker: ' + str(num)
  subprocess.Popen(['ssh', '-q', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + worker_ip, worker_command])


def stop_experiment(controller_ip, worker_ips):
  stop_controller(controller_ip);

  assert(config.WORKER_NUM <= len(worker_ips))
  num = 0;
  for ip in worker_ips:
    num += 1
    stop_worker(ip, num);


def stop_controller(controller_ip):
  controller_command =  'cd ' + NIMBUS_ROOT + ';'
  controller_command += 'scripts/stop-controller.sh &> /dev/null'

  # print '** Stopping controller: ' + controller_ip
  subprocess.Popen(['ssh', '-q', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip, controller_command])


def stop_worker(worker_ip, num):
  worker_command =  'cd ' + NIMBUS_ROOT + ';'
  worker_command += 'scripts/stop-workers.sh &> /dev/null'

  # print '** Stopping worker: ' + str(num)
  subprocess.Popen(['ssh', '-q', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + worker_ip, worker_command])


def test_nodes(node_ips):
  command =  'cd ' + NIMBUS_ROOT + ';'
  command += 'pwd;'

  num = 0;
  for ip in node_ips:
    num += 1;
    # print '** Testing node: ' + str(num) + ' ip: ' + ip
    subprocess.Popen(['ssh', '-q', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip, command])


def collect_logs(controller_ip, worker_ips):

  subprocess.call(['rm', '-rf', OUTPUT_PATH])
  subprocess.call(['mkdir', '-p', OUTPUT_PATH])

  collect_controller_logs(controller_ip)

  assert(config.WORKER_NUM <= len(worker_ips))
  num = 0;
  for ip in worker_ips:
    num += 1
    collect_worker_logs(ip, num);


def collect_controller_logs(controller_ip):

  subprocess.Popen(['scp', '-q', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip + ':' + NIMBUS_ROOT + STD_OUT_LOG,
      OUTPUT_PATH])

  subprocess.Popen(['scp', '-q', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip + ':' + NIMBUS_ROOT + REL_LOGS_PATH,
      OUTPUT_PATH])

  subprocess.Popen(['scp', '-q', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip + ':' + NIMBUS_ROOT +
      REL_CONTROLLER_PATH + LOAD_BALANCER_LOG,
      OUTPUT_PATH])

  subprocess.Popen(['scp', '-q', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip + ':' + NIMBUS_ROOT +
      REL_CONTROLLER_PATH + CONTROLLER_PER_ITER_STAT_LOG,
      OUTPUT_PATH])


def collect_worker_logs(worker_ip, num):

  subprocess.Popen(['scp', '-q', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + worker_ip + ':' + NIMBUS_ROOT + '*_' + STD_OUT_LOG,
      OUTPUT_PATH])

  subprocess.Popen(['scp', '-q', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + worker_ip + ':' + NIMBUS_ROOT + REL_LOGS_PATH,
      OUTPUT_PATH])

  subprocess.Popen(['scp', '-q', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + worker_ip + ':' + NIMBUS_ROOT +
      REL_WORKER_PATH + '*_' + WORKER_FINAL_STAT_LOG,
      OUTPUT_PATH])

  subprocess.Popen(['scp', '-q', '-r', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + worker_ip + ':' + NIMBUS_ROOT +
      REL_WORKER_PATH + '*_' + WORKER_PER_ITER_STAT_LOG,
      OUTPUT_PATH])

#  subprocess.Popen(['scp', '-q', '-r', '-i', config.PRIVATE_KEY,
#      '-o', 'UserKnownHostsFile=/dev/null',
#      '-o', 'StrictHostKeyChecking=no',
#      'ubuntu@' + worker_ip + ':' + NIMBUS_ROOT +
#      REL_WORKER_PATH + 'core',
#      OUTPUT_PATH + str(num) + '_core'])


def clean_logs(controller_ip, worker_ips):
  worker_path     = NIMBUS_ROOT + REL_WORKER_PATH;
  controller_path = NIMBUS_ROOT + REL_CONTROLLER_PATH;

  command  =  ''
  command +=  'rm -rf ' + NIMBUS_ROOT + REL_LOGS_PATH + ';'
  command +=  'rm -rf ' + NIMBUS_ROOT + '*' + STD_OUT_LOG + ';'

  command +=  'rm -rf ' + controller_path + '*.txt;'
  command +=  'rm -rf ' + controller_path + '*log*;'
  command +=  'rm -rf ' + controller_path + 'core;'

  command +=  'rm -rf ' + worker_path + '*.txt;'
  command +=  'rm -rf ' + worker_path + '*log*;'
  command +=  'rm -rf ' + worker_path + 'core;'
  command +=  'rm -rf ' + worker_path + '_db*;'
  command +=  'rm -rf ' + worker_path + 'split_output/;'
  command +=  'rm -rf ' + worker_path + 'output/;'

  # print '** Cleaning controller: ' + controller_ip
  subprocess.Popen(['ssh', '-q', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip, command])


  num = 0
  for ip in worker_ips:
    num += 1
  
    # print '** Cleaning worker ' + str(num)
    subprocess.Popen(['ssh', '-q',  '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip, command])


