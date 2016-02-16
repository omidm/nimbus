#!/usr/bin/python

import boto.ec2
import time
import sys

kRegionName = 'us-west-2'
kKeyName = 'quhang_temp_test'
kSecurityGroup = 'nimbus_sg_uswest2'
timestamp = sys.argv[1]
connection = boto.ec2.connect_to_region(kRegionName)
assert connection, 'Connection failed.'

instances = [i for i in connection.get_only_instances()
             if i.placement_group == timestamp]

while True:
    flag = True
    for instance in instances:
        print instance.ip_address, instance.instance_type, instance.state
        if instance.state != 'running':
            flag = False
    if flag:
        break
    print '-----------------------'
    time.sleep(5)

f = open(timestamp,'w')
for instance_type in ['c3.xlarge', 'c3.2xlarge', 'c3.4xlarge', 'c3.8xlarge']:
    ip_address = []
    for instance in instances:
        if str(instance.instance_type) == instance_type:
            ip_address.append(instance.ip_address)
    f.write(instance_type)
    f.write('\n')
    f.write(str(ip_address))
    f.write('\n')
    for i in ip_address:
        f.write('{} slots=8\n'.format(i))
f.close()

line = raw_input('Terminate yes?\n')
if line == 'yes':
    connection.delete_placement_group(timestamp)
    for instance in instances:
        instance.terminate()
