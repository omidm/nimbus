#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>


# US West (Northern California) Region
# EC2_LOCATION = 'us-west-1'
# NIMBUS_AMI = 'ami-50201815'
# UBUNTU_AMI = 'ami-660c3023'
# KEY_NAME = 'omidm-sing-key-pair-us-west-1'
# SG_GROUP = 'nimbus_sg_uswest1'
# INSTANCE_TYPE = 't1.micro'
# PLACEMENT_GROUP = None
# INSTANCE_NUM = 8


EC2_LOCATION = 'us-west-2'
NIMBUS_AMI = 'ami-5a45306a'
UBUNTU_AMI = 'ami-fa9cf1ca'
KEY_NAME = 'omidm-sing-key-pair-us-west-2'
SG_GROUP = 'nimbus_sg_uswest2'
# INSTANCE_TYPE = 't1.micro'
INSTANCE_TYPE = 'c3.xlarge'
PLACEMENT_GROUP = 'nimbus-cluster'

INSTANCE_NUM = 8


