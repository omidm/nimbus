#!/usr/bin/env python

import sys
import os
import re
import numpy
from decimal import *


sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import config



force_mean_time = []
force_std_time = []
advect_phi_mean_time = []
advect_phi_std_time = []
advect_removed_particles_mean_time = []
advect_removed_particles_std_time = []

for i in range(1, config.INSTANCE_NUM + 1):
  print 'Reading the file for sub grid ' + str(i)

  file_name = sys.argv[1] + '/' + str(i) + '/common/log.txt'
  f = open(file_name, 'r')
  content = f.readlines()
  
  force = []
  advect_phi = []
  advect_removed_particles = []
  
  collect = False
  for line in content:
    if not collect:
      result = re.findall('.*Calculate Dt.*', line)
      if len(result) > 0:
        force_time = 0
        advect_phi_time = 0
        advect_removed_particles_time = 0
        collect = True
        continue
    else:
      result = re.findall('.*Project.*', line)
      if len(result) > 0:
        force.append(force_time)
        advect_phi.append(advect_phi_time)
        advect_removed_particles.append(advect_removed_particles_time)
        collect = False
        continue
  
      result = re.findall('.*Forces\".*(\d+\.\d+)', line)
      if len(result) > 0:
        force_time += float(result[0]) * 1000
        continue
  
      result = re.findall('.*Advect Phi\".*(\d+\.\d+)', line)
      if len(result) > 0:
        advect_phi_time += float(result[0]) * 1000
        continue
  
      result = re.findall('.*Advect Removed Particles\".*(\d+\.\d+)', line)
      if len(result) > 0:
        advect_removed_particles_time += float(result[0]) * 1000
        continue
  
  f.close()
  
  
  force_mean_time.append(numpy.mean(force))
  force_std_time.append(numpy.std(force))
  
  advect_phi_mean_time.append(numpy.mean(advect_phi))
  advect_phi_std_time.append(numpy.std(advect_phi))
  
  advect_removed_particles_mean_time.append(numpy.mean(advect_removed_particles))
  advect_removed_particles_std_time.append(numpy.std(advect_removed_particles))






import numpy as np
import matplotlib.pyplot as plt
from decimal import *
getcontext().prec = 2

n_groups = config.INSTANCE_NUM
y_ticks = 3
index = np.arange(n_groups)

# mean_time = (np.mean(force), np.mean(advect_phi), np.mean(advect_removed_particles))
# std_time = (np.std(force), np.std(advect_phi), np.std(advect_removed_particles))

mean_time = force_mean_time
std_time = force_std_time

fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True)

error_config = {'ecolor': '1'}

ax0.errorbar(index, force_mean_time, yerr=force_std_time, fmt='o')
ax0.set_title('Apply Force Job')
# ax0.set_yscale('log')
ax0.locator_params(axis = 'y', nbins = y_ticks)

# y_min = min(force_mean_time) - max(force_std_time)
# y_max = max(force_mean_time) + max(force_std_time)
# y_step = (y_max - y_min) / float(y_ticks)
# ax0.yaxis.set_ticks(np.arange(y_min, y_max, y_step))

ax1.errorbar(index, advect_phi_mean_time, yerr=advect_phi_std_time, fmt='o')
ax1.set_title('Advect Phi Job')
# ax1.set_yscale('log')
ax1.locator_params(axis = 'y', nbins = y_ticks)

ax2.errorbar(index, advect_removed_particles_mean_time, yerr=advect_removed_particles_std_time, fmt='o')
ax2.set_title('Advect Removed Particles Job')
# ax2.set_yscale('log')
ax2.locator_params(axis = 'y', nbins = y_ticks)

plt.xlim((-1, n_groups))
plt.xlabel('Partition')
plt.xticks(index, ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII'))

# plt.grid()
plt.savefig("test.png")
plt.show()












def find_version_table(file_name, tag):
  f = open(file_name, 'r')
  content = f.readlines();
 
  count = 0;
  table = {}
  regexp =  '.*' + tag + '.*Compute:(\w+)\s*id:\s*(\d+)\s*version_hash:\s*(\d+).*'
  for line in content:
    result = re.findall(regexp, line)
    if (len(result) >= 1):
      table[(result[0][0] + "-" + result[0][1])] = result[0][2]
      count = count + 1
    # else:
    #   print "ERROR: corrupted line in " + file_name
    #   print line
  print "Read " + str(count) + " valid lines in " + file_name
  return table


def compare_version_tables(table_ref, table):
  match = True
  for key in table.keys():
    if not table_ref.has_key(key):
      print "ERROR: could not find key in refernce table: " + key
      match = False
    elif table_ref[key] != table[key]:
      print "ERROR: the value does not match for the key: " + key
      match = False
  return match




