#!/bin/bash

mkdir -p ~/.ssh 
touch ~/.ssh/authorized_keys
cat ~/cloud/src/nimbus/scripts/nimbus-rsa-key-pair.pub >> ~/.ssh/authorized_keys

