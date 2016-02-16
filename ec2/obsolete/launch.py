#!/usr/bin/python

import boto.ec2
import time
import sys

def filterImage(connection, image_name):
    images = [image for image in connection.get_all_images(
            owners='self',filters={'name':image_name})]
    image_names = [image.name for image in images]
    return (images, image_names)

kRegionName = 'us-west-2'
kKeyName = 'quhang_temp_test'
kSecurityGroup = 'nimbus_sg_uswest2'
timestamp = ''.join('_' if c in {' ', ':', '/'} else c
                    for c in str(time.strftime('%X %x %Z')))
connection = boto.ec2.connect_to_region(kRegionName)
assert connection, 'Connection failed.'

history = dict()
while True:
    line = raw_input()
    if line == '.':
        break
    images, image_names = filterImage(connection, line)
    print image_names
assert len(images) == 1, 'Not valid image.'
history['image_id'] = images[0]
history['image_name'] = image_names[0]
print history

stat = connection.create_placement_group(timestamp)
assert stat, 'Cannot build placement group'

instances = []
while True:
    line = raw_input('x/2x/4x/8x [num of machine]\n')
    if line == '.':
        break
    items = line.split()
    instance_type = 'c3.{}large'.format(items[0])
    instance_num = int(items[1])
    print instance_type, '#machine', instance_num, 'Confirm?'
    line = raw_input()
    if line == 'yes':
        reservation = history['image_id'].run(
                min_count=instance_num,
                max_count=instance_num,
                key_name=kKeyName,
                security_groups=[kSecurityGroup],
                placement_group=timestamp,
                instance_type=instance_type)
        instances.extend(reservation.instances)
        continue
