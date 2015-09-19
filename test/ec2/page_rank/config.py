#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

# ssh -i ~/.ssh/omidm-sing-key-pair-us-west-2.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ubuntu@<ip>

# US West (Northern California) Region
# EC2_LOCATION = 'us-west-1'
# NIMBUS_AMI = 'ami-50201815'
# UBUNTU_AMI = 'ami-660c3023'
# KEY_NAME = 'omidm-sing-key-pair-us-west-1'
# SECURITY_GROUP = 'nimbus_sg_uswest1'
# # INSTANCE_TYPE = 't1.micro'
# INSTANCE_TYPE = 'c3.xlarge'
# PLACEMENT_GROUP = None
# INSTANCE_NUM = 3

# EC2 configurations
# US West (Oregon) Region
EC2_LOCATION                    = 'us-west-2'
UBUNTU_AMI                      = 'ami-fa9cf1ca'
NIMBUS_AMI                      = 'ami-37202e07'
KEY_NAME                        = 'nimbus_chinmayee'
SECURITY_GROUP                  = 'nimbus_sg_uswest2'
CONTROLLER_INSTANCE_TYPE        = 'c3.4xlarge'
WORKER_INSTANCE_TYPE            = 'c3.2xlarge'
PLACEMENT                       = 'us-west-2c' # None
PRIVATE_KEY                     = '/Users/chinmayee/.ssh/' + KEY_NAME + '.pem'
CONTROLLER_NUM                  = 1
WORKER_NUM                      = 5
INSTANCE_NUM                    = WORKER_NUM + CONTROLLER_NUM

# controller configurations
ASSIGNER_THREAD_NUM             = 8
BATCH_ASSIGN_NUM                = 200
COMMAND_BATCH_SIZE              = 10000
DEACTIVATE_CONTROLLER_TEMPLATE  = False  # NOTE: CHANGE THIS FOR EXPERIMENTS
DEACTIVATE_COMPLEX_MEMOIZATION  = False
DEACTIVATE_BINDING_MEMOIZATION  = False
DEACTIVATE_WORKER_TEMPLATE      = False
DEACTIVATE_MEGA_RCR_JOB         = False
DEACTIVATE_DM_QUERY_CACHE       = False
ACTIVATE_LB                     = False
ACTIVATE_FT                     = False
LB_PERIOD                       = 60
FT_PERIOD                       = 600
FIRST_PORT                      = 5900


# worker/simulation (app specific) configurations
PARTITIONS                      = 80  # NOTE: CHANGE THIS FOR EXPERIMENTS
INPUT_DIR                       = 'input-wiki-' + str(PARTITIONS) + 'p'
OUTPUT_DIR                      = 'output-wiki'
ITERATION_NUM                   = 11


# placement group
# PLACEMENT_GROUP             = '*'
PLACEMENT_GROUP               = 'pagerank-' + str(PARTITIONS) + 'p-' + str(WORKER_NUM) + 'w'  # None


# worker/simulation (common) configurations
# OTHREAD_NUM                     = (PARTITIONS/WORKER_NUM)
OTHREAD_NUM                     = 8  # NOTE: VERIFY THAT THIS MAKES SENSE
RUN_WITH_TASKSET                = False
WORKER_TASKSET                  = '0-1,4-5'
# ITERATION_BATCH               = 1


# logging configurations
STD_OUT_LOG                     = 'ec2_log.txt'
LOAD_BALANCER_LOG               = 'load_balancer_log'
SCHED_PER_ITER_STAT_LOG         = 'controller_stats.txt'
WORKER_FINAL_STAT_LOG           = 'time_per_thread.txt' 
WORKER_PER_ITER_STAT_LOG        = 'main_timers.txt' 
OUTPUT_PATH                     = 'output/'


# Path configuration
SOURCE_NIMBUS_ROOT   = '~/cloud/src/nimbus/'
EC2_NIMBUS_ROOT      = '~/cloud/src/nimbus/'
REL_APPLICATION_PATH = 'application/page_rank'
REL_SCHEDULER_PATH   = 'test/scheduler_v3/'
REL_WORKER_PATH      = 'test/page_rank/'

