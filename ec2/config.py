#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

# ssh -i ~/.ssh/omidm-sing-key-pair-us-west-2.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ubuntu@<ip>

# US West (Northern California) Region
# EC2_LOCATION = 'us-west-1'
# NIMBUS_AMI = 'ami-50201815'
# UBUNTU_AMI = 'ami-660c3023'
# KEY_NAME = 'omidm-sing-key-pair-us-west-1'
# SECURITY_GROUP = 'nimbus_sg_uswest1'

# EC2 configurations
# US West (Oregon) Region
EC2_LOCATION                    = 'us-west-2'
UBUNTU_AMI                      = 'ami-fa9cf1ca'
NIMBUS_AMI                      = 'ami-66af4c06'
CONTROLLER_INSTANCE_TYPE        = 'c3.4xlarge'
WORKER_INSTANCE_TYPE            = 'c3.2xlarge'
PLACEMENT                       = 'us-west-2c' # None
PLACEMENT_GROUP                 = 'nimbus-cluster' # None / '*'
SECURITY_GROUP                  = 'nimbus_sg_uswest2'
KEY_NAME                        = 'omidm-sing-key-pair-us-west-2'
PRIVATE_KEY                     = '/home/omidm/.ssh/' + KEY_NAME + '.pem'
CONTROLLER_NUM                  = 1
WORKER_NUM                      = 2


#Environment variables
DBG_MODE                        = 'error'
TTIMER_LEVEL                    = 'l1' 


# Controller configurations
ASSIGNER_THREAD_NUM             = 8
BATCH_ASSIGN_NUM                = 200
COMMAND_BATCH_SIZE              = 10000
DEACTIVATE_CONTROLLER_TEMPLATE  = False
DEACTIVATE_COMPLEX_MEMOIZATION  = False
DEACTIVATE_BINDING_MEMOIZATION  = False
DEACTIVATE_WORKER_TEMPLATE      = False
DEACTIVATE_MEGA_RCR_JOB         = False
DEACTIVATE_CASCADED_BINDING     = False
ACTIVATE_LB                     = False
ACTIVATE_FT                     = False
LB_PERIOD                       = 60
FT_PERIOD                       = 600
FIRST_PORT                      = 5800
SPLIT_ARGS                      = str(WORKER_NUM) + ' 1 1'


# Worker configurations
OTHREAD_NUM                     = 8
APPLICATION                     = 'lr' # 'k-means' 'water' 
RUN_WITH_TASKSET                = False
WORKER_TASKSET                  = '0-1,4-5' # '0-3,8-11'


# Application configurations
DIMENSION                       = 10
ITERATION_NUM                   = 30
PARTITION_NUM                   = 2000
SAMPLE_NUM_M                    = 100


# LR
LR_APP_OPTIONS  = ' '
LR_APP_OPTIONS += ' -d ' + str(DIMENSION)
LR_APP_OPTIONS += ' -i ' + str(ITERATION_NUM)
LR_APP_OPTIONS += ' --sn ' + str(SAMPLE_NUM_M)
LR_APP_OPTIONS += ' --pn ' + str(PARTITION_NUM)
LR_REL_APP_PATH = 'applications/ml/logistic_regression/liblogistic_regression.so' 



