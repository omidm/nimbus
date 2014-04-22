#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>


# US West (Northern California) Region
# EC2_LOCATION = 'us-west-1'
# NIMBUS_AMI = 'ami-50201815'
# UBUNTU_AMI = 'ami-660c3023'
# KEY_NAME = 'omidm-sing-key-pair-us-west-1'
# SECURITY_GROUP = 'nimbus_sg_uswest1'
# INSTANCE_TYPE = 't1.micro'
# PLACEMENT_GROUP = None
# INSTANCE_NUM = 1


# US West (Oregon) Region
EC2_LOCATION = 'us-west-2'
NIMBUS_AMI = 'ami-5a45306a'
UBUNTU_AMI = 'ami-fa9cf1ca'
KEY_NAME = 'omidm-sing-key-pair-us-west-2'
SECURITY_GROUP = 'nimbus_sg_uswest2'
# INSTANCE_TYPE = 't1.micro'
# PLACEMENT_GROUP = None
INSTANCE_TYPE = 'c3.xlarge'
PLACEMENT_GROUP = 'nimbus-cluster'
INSTANCE_NUM = 8




PRIVATE_KEY = '/home/omidm/.ssh/' + KEY_NAME + '.pem'
DIRECTORY_PATH = 'test/src/physbam/'
NODES_FILE_NAME = 'test_nodes.txt' 
OUTPUT_NAME = 'output/' 
SCALE = 256
FRAME_NUM = 40
