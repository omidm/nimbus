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
NIMBUS_ROOT                     = '~/cloud/src/nimbus/'
REL_LOGS_PATH                    = 'logs/'
REL_CONTROLLER_PATH             = 'nodes/nimbus_controller/'
REL_WORKER_PATH                 = 'nodes/nimbus_worker/'


def start_experiment(controller_ip, controller_p_ip, worker_ips, worker_p_ips):
  worker_num = len(worker_ips)
  start_controller(controller_ip, worker_num);
  time.sleep(5)

  idx = 0;
  for idx in range(0, len(worker_ips)):
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

  print '** Starting controller: ' + controller_ip
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
  if config.RUN_WITH_TASKSET:
    worker_command += 'taskset -c ' + config.WORKER_TASKSET + './scripts/start-workers.sh 1'
  else:
    worker_command += './scripts/start-workers.sh 1'

  worker_command += ' -p ' + str(config.FIRST_PORT + num)
  worker_command += ' --ip ' + worker_p_ip
  worker_command += ' --cip ' + controller_p_ip
  worker_command += ' --cport ' + str(config.FIRST_PORT)
  worker_command += ' --othread ' + str(config.OTHREAD_NUM)
  if (config.APPLICATION == 'lr'):
    worker_command += ' -l ' + str(config.LR_REL_APP_PATH)
    worker_command += ' ' + str(config.LR_APP_OPTIONS)
  else:
    print "ERROR: unknown application: " + config.APPLICATION
    exit(0)
  worker_command += ' &>> ' + str(num) + '_' + STD_OUT_LOG

  print '** Starting worker: ' + str(num)
  subprocess.Popen(['ssh', '-q', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + worker_ip, worker_command])


def stop_experiment(controller_ip, worker_ips):
  stop_controller(controller_ip);

  num = 0;
  for ip in worker_ips:
    num += 1
    stop_worker(ip, num);


def stop_controller(controller_ip):
  controller_command =  'cd ' + NIMBUS_ROOT + ';'
  controller_command += 'scripts/stop-controller.sh &> /dev/null'

  print '** Stopping controller: ' + controller_ip
  subprocess.Popen(['ssh', '-q', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip, controller_command])


def stop_worker(worker_ip, num):
  worker_command =  'cd ' + NIMBUS_ROOT + ';'
  worker_command += 'scripts/stop-workers.sh &> /dev/null'

  print '** Stopping worker: ' + str(num)
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
    print '** Testing node: ' + str(num) + ' ip: ' + ip
    subprocess.Popen(['ssh', '-q', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip, command])


def collect_logs(controller_ip, worker_ips):

  subprocess.call(['rm', '-rf', OUTPUT_PATH])
  subprocess.call(['mkdir', '-p', OUTPUT_PATH])

  collect_controller_logs(controller_ip)

  num = 0
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

  subprocess.Popen(['ssh', '-q', '-i', config.PRIVATE_KEY,
      '-o', 'UserKnownHostsFile=/dev/null',
      '-o', 'StrictHostKeyChecking=no',
      'ubuntu@' + controller_ip, command])

  print '** Cleaning controller ' + controller_ip

  num = 0
  for ip in worker_ips:
    num += 1
  
    print '** Cleaning worker ' + str(num)
    subprocess.Popen(['ssh', '-q',  '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip, command])


