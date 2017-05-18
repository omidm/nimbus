#!/bin/bash

# nimbus
sudo apt-get update
sudo apt-get install git
sudo apt-get install gcc
sudo apt-get install g++
sudo apt-get install make
sudo apt-get install python

# for plotting scripts
sudo apt-get install python-numpy
sudo apt-get install python-matplotlib

# for visualizing physbam through scripts
sudo apt-get install xdotool

# for compiling the mark down docs
sudo apt-get install pandoc

# for creating fake stragglers with cgroup
sudo apt-get install cgroup-bin

# to use the ec2 scripts
sudo apt install python-pip
pip install -U boto
