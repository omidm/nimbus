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
  time.sleep(5)

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
  scheduler_command += ' --split ' + str(config.PARTITIONS) + ' 1 1'
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

  # app specific
  worker_command += ' --input ' + str(config.INPUT_DIR)
  worker_command += ' --output ' + str(config.OUTPUT_DIR)
  worker_command += ' --iterations ' + str(config.ITERATION_NUM)

  # general
  worker_command += ' --othread ' + str(config.OTHREAD_NUM)
  # worker_command += ' --ibatch ' + str(config.ITERATION_BATCH)
  worker_command += ' --port ' + str(config.FIRST_PORT + num)
  worker_command += ' --ip ' + worker_p_ip
  worker_command += ' --cip ' + scheduler_p_ip
  worker_command += ' --cport ' + str(config.FIRST_PORT)
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
        config.REL_WORKER_PATH + '*_' + config.WORKER_FINAL_STAT_LOG,
        config.OUTPUT_PATH])

    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
        '-o', 'UserKnownHostsFile=/dev/null',
        '-o', 'StrictHostKeyChecking=no',
        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
        config.REL_WORKER_PATH + '*_' + config.WORKER_PER_ITER_STAT_LOG,
        config.OUTPUT_PATH])

#    subprocess.Popen(['scp', '-r', '-i', config.PRIVATE_KEY,
#        '-o', 'UserKnownHostsFile=/dev/null',
#        '-o', 'StrictHostKeyChecking=no',
#        'ubuntu@' + ip + ':' + config.EC2_NIMBUS_ROOT +
#        config.REL_WORKER_PATH + 'core',
#        config.OUTPUT_PATH + str(num) + '_core'])


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


