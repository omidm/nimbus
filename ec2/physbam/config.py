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
# INSTANCE_NUM = 4


# US West (Oregon) Region
EC2_LOCATION = 'us-west-2'
NIMBUS_AMI = 'ami-5a45306a'
UBUNTU_AMI = 'ami-fa9cf1ca'
KEY_NAME = 'omidm-sing-key-pair-us-west-2'
SECURITY_GROUP = 'nimbus_sg_uswest2'
INSTANCE_TYPE = 'c3.2xlarge'
PLACEMENT_GROUP = 'physbam-cluster'
INSTANCE_NUM = 8


FOLDER_PATH = '~/physbam/Projects/Water/'
OUTPUT_PATH = "output/"


PRIVATE_KEY = '/home/omidm/.ssh/' + KEY_NAME + '.pem'
REMOTE_PATH = '~/test/src/physbam/'
SOURCE_PATH = '~/cloud/src/nimbus/application/physbam-app/water/'
NODES_FILE_NAME = 'test_nodes.txt' 
OUTPUT_NAME = 'output/' 
SCALE = 256
FRAME_NUM = 40
