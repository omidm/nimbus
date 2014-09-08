#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>


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
EC2_LOCATION     = 'us-west-2'
UBUNTU_AMI       = 'ami-fa9cf1ca'
NIMBUS_0_3_AMI   = 'ami-69ee9559'
KEY_NAME         = 'omidm-sing-key-pair-us-west-2'
SECURITY_GROUP   = 'nimbus_sg_uswest2'
INSTANCE_TYPE    = 'c3.2xlarge' # 't1.micro'
PLACEMENT_GROUP  = 'nimbus-cluster' # None
PRIVATE_KEY      = '/home/omidm/.ssh/' + KEY_NAME + '.pem'

# Experiment configurations
WORKER_NUM           = 100
SCHEDULER_NUM        = 1
ASSIGNER_THREAD_NUM  = 4
BATCH_ASSIGN_NUM     = 200
OTHREAD_NUM          = 4
INSTANCE_NUM         = WORKER_NUM + SCHEDULER_NUM
FIRST_PORT           = 5800
LOG_FILE_NAME        = 'ec2_log.txt'
SCHED_LOG_FILE_NAME  = 'job_assigner_log'
WORKER_LOG_FILE_NAME = 'worker-log-'
OUTPUT_PATH          = 'output/'


# Build and Run configuration
SOURCE_NIMBUS_ROOT   = '~/cloud/src/nimbus/'
EC2_NIMBUS_ROOT      = '~/cloud/src/nimbus/'
REL_APPLICATION_PATH = 'application/water_multiple/Build/Release/'
REL_SCHEDULER_PATH   = 'test/scheduler_v2/'
REL_WORKER_PATH      = 'test/water_multiple/'

NIMBUS_LIB           = 'libnimbus.so'
APPLICATION_LIB      = 'libwater_app.so'
SCHEDULER_BINARY     = 'scheduler'
WORKER_BINARY        = 'worker'

