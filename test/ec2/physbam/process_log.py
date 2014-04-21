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
from matplotlib.font_manager import FontProperties
from decimal import *
getcontext().prec = 2

n_groups = config.INSTANCE_NUM
index = np.arange(n_groups)
y_ticks = 3
font = FontProperties()
font.set_family('serif')

fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True)

for i in np.arange(len(force_std_time)):
  if (force_mean_time[i] <= force_std_time[i]):
    force_std_time[i] = force_mean_time[i] - 1

ax0.errorbar(index, force_mean_time, yerr=force_std_time, fmt='o')
ax0.set_title('Apply Force Job', fontproperties=font)
# ax0.set_yscale('log')
# ax0.set_ylabel('log time(ms)', fontproperties=font)
ax0.locator_params(axis = 'y', nbins = y_ticks)
ax0.set_ylabel('time(ms)', fontproperties=font)


for i in np.arange(len(advect_phi_std_time)):
  if (advect_phi_mean_time[i] <= advect_phi_std_time[i]):
    advect_phi_std_time[i] = advect_phi_mean_time[i] - 1

ax1.errorbar(index, advect_phi_mean_time, yerr=advect_phi_std_time, fmt='o')
ax1.set_title('Advect Phi Job', fontproperties=font)
# ax1.set_yscale('log')
# ax1.set_ylabel('log time(ms)', fontproperties=font)
ax1.locator_params(axis = 'y', nbins = y_ticks)
ax1.set_ylabel('time(ms)', fontproperties=font)


for i in np.arange(len(advect_removed_particles_std_time)):
  if (advect_removed_particles_mean_time[i] <= advect_removed_particles_std_time[i]):
    advect_removed_particles_std_time[i] = advect_removed_particles_mean_time[i] - 1

ax2.errorbar(index, advect_removed_particles_mean_time, yerr=advect_removed_particles_std_time, fmt='o')
ax2.set_title('Advect Removed Particles Job', fontproperties=font)
ax2.set_yscale('log')
ax2.set_ylabel('log time(ms)', fontproperties=font)
# ax2.locator_params(axis = 'y', nbins = y_ticks)
# ax2.set_ylabel('time(ms)', fontproperties=font)

plt.xlim((-1, n_groups))
plt.xlabel('Partition', fontproperties=font)
plt.xticks(index, ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII'), fontproperties=font)

# plt.grid()
plt.savefig("test.png")
plt.show()









