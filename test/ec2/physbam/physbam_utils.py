#!/usr/bin/env python

# Author: Omid Mashayekhi <omidm@stanford.edu>

import sys
import os
import subprocess

import physbam_config

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import config

temp_file_name = '_temp_file_'

def make_nodes_file_content(ip_addresses):

  string = ""
  for ip in ip_addresses:
    string = string + ip + " cpu=1\n"

  file = open(temp_file_name, 'w+')
  file.write(string)
  file.close()



def copy_nodes_file_to_hosts(ip_addresses):

  make_nodes_file_content(ip_addresses)

  for ip in ip_addresses:
    subprocess.call(['scp', '-i', '/home/omidm/.ssh/' + config.KEY_NAME + '.pem',
        temp_file_name, 'ubuntu@' + ip + ':~/' + physbam_config.DIRECTORY_PATH +
        physbam_config.NODES_FILE_NAME])


def run_experiment(ip):
  print "omid"


  command =  'cd ' + physbam_config.DIRECTORY_PATH + ';'
  command += 'mpirun -hostfile ' + physbam_config.NODES_FILE_NAME
  command += ' -np ' + str(config.INSTANCE_NUM)
  command += ' ./Water -scale 20 -e 10'

  subprocess.call(['ssh', '-i', '/home/omidm/.ssh/' + config.KEY_NAME + '.pem',
      'ubuntu@' + ip, command])
