#!/usr/bin/env python

import sys
import os
import re
import numpy
from decimal import *

import config

force_agg = 0;
advect_phi_agg = 0;
advect_v_agg = 0;
advect_removed_particles_agg = 0;

exchange_force_agg = 0;
exchange_advect_phi_agg = 0;
exchange_advect_v_agg = 0;
exchange_advect_removed_particles_agg = 0;



force_mean_time = []
force_std_time = []
advect_phi_mean_time = []
advect_phi_std_time = []
advect_v_mean_time = []
advect_v_std_time = []
advect_removed_particles_mean_time = []
advect_removed_particles_std_time = []

for i in range(1, config.INSTANCE_NUM + 1):
  print 'Reading the file for sub grid ' + str(i)

  file_name = sys.argv[1] + '/' + str(i) + '/common/log.txt'
  f = open(file_name, 'r')
  content = f.readlines()
  
  force = []
  advect_phi = []
  advect_v = []
  advect_removed_particles = []
  
  collect = False
  for line in content:
    if not collect:
      result = re.findall('.*Calculate Dt.*', line)
      if len(result) > 0:
        force_time = 0
        advect_phi_time = 0
        advect_v_time = 0
        advect_removed_particles_time = 0
        collect = True
        continue
    else:
      result = re.findall('.*Project.*', line)
      if len(result) > 0:
        force.append(force_time)
        advect_phi.append(advect_phi_time)
        advect_v.append(advect_v_time)
        advect_removed_particles.append(advect_removed_particles_time)
        collect = False
        continue
  
      result = re.findall('.*Forces\".*(\d+\.\d+)', line)
      if len(result) > 0:
        force_time += float(result[0]) * 1000
        force_agg += float(result[0])
        continue
  
      result = re.findall('.*Advect Phi\".*(\d+\.\d+)', line)
      if len(result) > 0:
        advect_phi_time += float(result[0]) * 1000
        advect_phi_agg += float(result[0])
        continue
  
      result = re.findall('.*Advect V\".*(\d+\.\d+)', line)
      if len(result) > 0:
        advect_v_time += float(result[0]) * 1000
        advect_v_agg += float(result[0])
        continue
 
      result = re.findall('.*Advect Removed Particles\".*(\d+\.\d+)', line)
      if len(result) > 0:
        advect_removed_particles_time += float(result[0]) * 1000
        advect_removed_particles_agg += float(result[0])
        continue
   
      result = re.findall('.*Forces Exchange\".*(\d+\.\d+)', line)
      if len(result) > 0:
        exchange_force_agg += float(result[0])
        continue
  
      result = re.findall('.*Advect Phi Exchange\".*(\d+\.\d+)', line)
      if len(result) > 0:
        exchange_advect_phi_agg += float(result[0])
        continue
   
      result = re.findall('.*Advect V Exchange\".*(\d+\.\d+)', line)
      if len(result) > 0:
        exchange_advect_v_agg += float(result[0])
        continue
  
      result = re.findall('.*Advect Removed Particles Exchange\".*(\d+\.\d+)', line)
      if len(result) > 0:
        exchange_advect_removed_particles_agg += float(result[0])
        continue
 
  f.close()
  
  
  force_mean_time.append(numpy.mean(force))
  force_std_time.append(numpy.std(force))
  
  advect_phi_mean_time.append(numpy.mean(advect_phi))
  advect_phi_std_time.append(numpy.std(advect_phi))
  
  advect_v_mean_time.append(numpy.mean(advect_v))
  advect_v_std_time.append(numpy.std(advect_v))
  
  advect_removed_particles_mean_time.append(numpy.mean(advect_removed_particles))
  advect_removed_particles_std_time.append(numpy.std(advect_removed_particles))



for i in numpy.arange(len(force_std_time)):
  if (force_mean_time[i] <= force_std_time[i]):
    force_std_time[i] = force_mean_time[i] - 1


for i in numpy.arange(len(advect_phi_std_time)):
  if (advect_phi_mean_time[i] <= advect_phi_std_time[i]):
    advect_phi_std_time[i] = advect_phi_mean_time[i] - 1

for i in numpy.arange(len(advect_v_std_time)):
  if (advect_v_mean_time[i] <= advect_v_std_time[i]):
    advect_v_std_time[i] = advect_v_mean_time[i] - 1

for i in numpy.arange(len(advect_removed_particles_std_time)):
  if (advect_removed_particles_mean_time[i] <= advect_removed_particles_std_time[i]):
    advect_removed_particles_std_time[i] = advect_removed_particles_mean_time[i] - 1



string = ''
string += 'Start\n'
string += 'd: '
for data in force_mean_time:
  string += str(data) + ' '
string += '\ne: '
for error in force_std_time:
  string += str(error) + ' '
string += '\nt: Apply Force Job\n'
string += 'End\n'

string += 'Start\n'
string += 'd: '
for data in advect_phi_mean_time:
  string += str(data) + ' '
string += '\ne: '
for error in advect_phi_std_time:
  string += str(error) + ' '
string += '\nt: Advect Phi Job\n'
string += 'End\n'

string += 'Start\n'
string += 'd: '
for data in advect_v_mean_time:
  string += str(data) + ' '
string += '\ne: '
for error in advect_v_std_time:
  string += str(error) + ' '
string += '\nt: Advect V Job\n'
string += 'End\n'


string += 'Start\n'
string += 'd: '
for data in advect_removed_particles_mean_time:
  string += str(data) + ' '
string += '\ne: '
for error in advect_removed_particles_std_time:
  string += str(error) + ' '
string += '\nt: Advect Removed Particles Job\n'
string += 'End\n'


string += "force communication: " + str(exchange_force_agg / float(force_agg))
string += "advect phi communication: " + str(exchange_advect_phi_agg / float(advect_phi_agg))
string += "advect v communication: " + str(exchange_advect_v_agg / float(advect_v_agg))
string +=  "advect removed particles communication: " + str(exchange_advect_removed_particles_agg / float(advect_removed_particles_agg))


f = open(sys.argv[2], 'w+')
f.write(string)
f.close()


# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.font_manager import FontProperties
# from decimal import *
# getcontext().prec = 2
# 
# n_groups = config.INSTANCE_NUM
# index = np.arange(n_groups)
# y_ticks = 3
# font = FontProperties()
# font.set_family('serif')
# 
# fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True)
# 
# for i in np.arange(len(force_std_time)):
#   if (force_mean_time[i] <= force_std_time[i]):
#     force_std_time[i] = force_mean_time[i] - 1
# 
# ax0.errorbar(index, force_mean_time, yerr=force_std_time, fmt='o')
# ax0.set_title('Apply Force Job', fontproperties=font)
# # ax0.set_yscale('log')
# # ax0.set_ylabel('log time(ms)', fontproperties=font)
# ax0.locator_params(axis = 'y', nbins = y_ticks)
# ax0.set_ylabel('time(ms)', fontproperties=font)
# 
# 
# for i in np.arange(len(advect_phi_std_time)):
#   if (advect_phi_mean_time[i] <= advect_phi_std_time[i]):
#     advect_phi_std_time[i] = advect_phi_mean_time[i] - 1
# 
# print  advect_phi_mean_time
# 
# ax1.errorbar(index, advect_phi_mean_time, yerr=advect_phi_std_time, fmt='o')
# ax1.set_title('Advect Phi Job', fontproperties=font)
# # ax1.set_yscale('log')
# # ax1.set_ylabel('log time(ms)', fontproperties=font)
# ax1.locator_params(axis = 'y', nbins = y_ticks)
# ax1.set_ylabel('time(ms)', fontproperties=font)
# 
# 
# for i in np.arange(len(advect_removed_particles_std_time)):
#   if (advect_removed_particles_mean_time[i] <= advect_removed_particles_std_time[i]):
#     advect_removed_particles_std_time[i] = advect_removed_particles_mean_time[i] - 1
# 
# ax2.errorbar(index, advect_removed_particles_mean_time, yerr=advect_removed_particles_std_time, fmt='o')
# ax2.set_title('Advect Removed Particles Job', fontproperties=font)
# # ax2.set_yscale('log')
# # ax2.set_ylabel('log time(ms)', fontproperties=font)
# ax2.locator_params(axis = 'y', nbins = y_ticks)
# ax2.set_ylabel('time(ms)', fontproperties=font)
# 
# plt.xlim((-1, n_groups))
# plt.xlabel('Partition', fontproperties=font)
# plt.xticks(index, ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII'), fontproperties=font)
# 
# # plt.grid()
# plt.savefig("test.png")
# plt.show()









