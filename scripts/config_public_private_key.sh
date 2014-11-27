#!/bin/bash

mkdir -p ~/.ssh 
touch ~/.ssh/authorized_keys

# create the key value pair
ssh-keygen -t rsa -N "" -f nimbus-rsa-key-pair

cat ~/cloud/src/nimbus/scripts/nimbus-rsa-key-pair.pub >> ~/.ssh/authorized_keys

