#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>


# US West (Northern California) Region
EC2_LOCATION = 'us-west-1'
NIMBUS_AMI = 'ami-50201815'
UBUNTU_AMI = 'ami-660c3023'
KEY_NAME = 'omidm-sing-key-pair-us-west-1'
SECURITY_GROUP = 'nimbus_sg_uswest1'
# INSTANCE_TYPE = 't1.micro'
INSTANCE_TYPE = 'c3.xlarge'
PLACEMENT_GROUP = None
INSTANCE_NUM = 3


# US West (Oregon) Region
# EC2_LOCATION = 'us-west-2'
# NIMBUS_AMI = 'ami-5a45306a'
# UBUNTU_AMI = 'ami-fa9cf1ca'
# KEY_NAME = 'omidm-sing-key-pair-us-west-2'
# SECURITY_GROUP = 'nimbus_sg_uswest2'
# INSTANCE_TYPE = 't1.micro'
# PLACEMENT_GROUP = None
# # INSTANCE_TYPE = 'c3.xlarge'
# # PLACEMENT_GROUP = 'nimbus-cluster'
# 
# INSTANCE_NUM = 2

PRIVATE_KEY = '/home/omidm/.ssh/' + KEY_NAME + '.pem'

SOURCE_NIMBUS_ROOT = '/home/omidm/cloud/src/nimbus/'
REL_APPLICATION_PATH = 'application/stencil_1d/'
REL_SCHEDULER_PATH = 'test/scheduler_1d/'
REL_WORKER_PATH = 'test/stencil_1d/'

NIMBUS_LIB = 'libnimbus.so'
APPLICATION_LIB = 'libstencil_1d.so'
SCHEDULER_BINARY = 'scheduler'
WORKER_BINARY = 'worker'

EC2_FOLDER_NAME = 'nimbus/'

FIRST_PORT = 5900
LOG_FILE_NAME = 'ec2_log.txt'
OUTPUT_PATH = 'output/'



