#!/usr/bin/env python

import boto.ec2

import config
import utils


# s3 = boto.connect_s3()
# bucket = s3.create_bucket('test.omidm')  # bucket names must be unique
# key = bucket.new_key('examples/first_file.csv')
# key.set_contents_from_filename('xyz')
# key.set_acl('public-read')


ec2 = boto.ec2.connect_to_region(config.EC2_LOCATION)

print utils.get_dns_names()
ip = utils.get_ip_addresses()
print ip[0]
print ip[1]

utils.wait_for_instances_to_start()

instances = ec2.get_all_instances()

for inst in instances:
  print inst.id
  print " "



instances = ec2.get_only_instances()

for inst in instances:
  print inst.id
  print " "




reservations = ec2.get_all_reservations()

for r in reservations:
  print "*** Reservation:"
  print r.id
  instances = r.instances
  for inst in instances:
    print "*Instance:" 
    print inst.public_dns_name 
    print inst.instance_type
    print inst.placement 
    print inst.id
    print inst.state
    # if inst.id == 'i-32a09f6d':
    #   print "Stoping Instance: ..."
    #   ec2.stop_instances(instance_ids=[inst.id])
    # if inst.id == 'i-32a09f6d':
    #   print "Terminating Instance: ..."
    #   ec2.terminate_instances(instance_ids=[inst.id])

# Ubuntu 12.04
# ec2.run_instances(
#     'ami-660c3023',
#     key_name='omidm-sing-key-pair-us-west-1',
#     instance_type='t1.micro',
#     security_groups=['nimbus_sg_uswest1'])

# Nimbus
# ec2.run_instances(
#     'ami-50201815',
#     key_name='omidm-sing-key-pair-us-west-1',
#     instance_type='t1.micro',
#     security_groups=['nimbus_sg_uswest1'])


# security_groups = ec2.get_all_security_groups()
# 
# for sg in security_groups:
#   print "*** Security Group:"
#   print sg.name
#   print sg.rules




