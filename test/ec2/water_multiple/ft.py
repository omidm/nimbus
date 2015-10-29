#!/usr/bin/env python

import sys
import argparse
import os
import re
import numpy
import decimal
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from operator import add
from matplotlib import rc


rc('font', size=24)

## Parse the command line arguments ##
parser = argparse.ArgumentParser(description='Process log files.')
parser.add_argument(
    "-if", "--input_ft",
    dest="iffname",
    default="_socc_ft_wt/ec2_log.txt",
    help="ft input file name")
parser.add_argument(
    "-in", "--input_normal",
    dest="infname",
    default="_socc_nimbus_256_16_wot_8w/ec2_log.txt",
    help="normal input file name")


args = parser.parse_args()

ft_file_name =  args.iffname
normal_file_name = args.infname


print 'Opening the ft file ' + ft_file_name
#f = open(ft_file_name, 'r')
#ft_content = f.readlines()

ft_time = [];
ft_base_time = 0;


# How to count number od iterations from std::out of controller? 
# either count loop_iteration_part_two jobs: 
#     cat log | grep complex | grep "projection_loop_iteration_end\." -c
#                      +
#     cat omid | grep Picked | grep " loop_iteration_part_two\."
#
# or count loop_iteration jobs (loop_frame is not templetized):
#    cat omid | grep complex | grep "loop_iteration_part_two\." -c
#                     +
#    cat omid | grep Picked | grep " loop_iteration\." -c
#

# for line in ft_content:
#   result = re.findall('.*: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) complex .* loop_iteration_part_two\.', line)
#   if len(result) > 0:
#     ft_time.append(decimal.Decimal(result[0]))
#     continue
#   result = re.findall('.*: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) Picked .* loop_iteration\.', line)
#   if len(result) > 0:
#     ft_time.append(decimal.Decimal(result[0]))
#     continue
#   result = re.findall('.*: (\d+|\d+\.\d+|\d+e-\d+|\d+\.\d+e-\d+) Picked .* main\.', line)
#   if len(result) > 0:
#     ft_base_time = decimal.Decimal(result[0])
# 
# f.close()

for i in range (0, len(ft_time)):
  ft_time[i] = ft_time[i] - ft_base_time

ft_time = [decimal.Decimal('14.294813633'), decimal.Decimal('27.324031115'),
           decimal.Decimal('40.443511248'), decimal.Decimal('53.069815159'),
           decimal.Decimal('65.462562800'), decimal.Decimal('77.733314038'),
           decimal.Decimal('90.061267138'), decimal.Decimal('102.221193075'),
           decimal.Decimal('114.438697338'), decimal.Decimal('126.760795355'),
           decimal.Decimal('138.974315167'), decimal.Decimal('151.146574259'),
           decimal.Decimal('163.287139178'), decimal.Decimal('175.477228403'),
           decimal.Decimal('194.298757553'), decimal.Decimal('206.415055037'),
           decimal.Decimal('218.860237837'), decimal.Decimal('231.237745285'),
           decimal.Decimal('243.899189711'), decimal.Decimal('256.256145716'),
           decimal.Decimal('268.684834242'), decimal.Decimal('281.248857975'),
           decimal.Decimal('293.468629837'), decimal.Decimal('305.821554661'),
           decimal.Decimal('318.039976120'), decimal.Decimal('330.304618359'),
           decimal.Decimal('342.729818583'), decimal.Decimal('355.131979943'),
           decimal.Decimal('367.641402483'), decimal.Decimal('380.125743389'),
           decimal.Decimal('392.634788990'), decimal.Decimal('405.396464825'),
           decimal.Decimal('424.999366761'), decimal.Decimal('437.769950867'),
           decimal.Decimal('450.453829050'), decimal.Decimal('463.325127602'),
           decimal.Decimal('475.998623371'), decimal.Decimal('488.751348973'),
           decimal.Decimal('501.526008368'), decimal.Decimal('514.184755087'),
           decimal.Decimal('526.848981619'), decimal.Decimal('539.389060498'),
           decimal.Decimal('552.218667031'), decimal.Decimal('564.810718775'),
           decimal.Decimal('577.659511328'), decimal.Decimal('609.571983815'),
           decimal.Decimal('623.033124924'), decimal.Decimal('635.865169525'),
           decimal.Decimal('648.381989718'), decimal.Decimal('661.265377760'),
           decimal.Decimal('674.238375664'), decimal.Decimal('686.919780016'),
           decimal.Decimal('841.143443346'), decimal.Decimal('866.256693125'),
           decimal.Decimal('884.582519055'), decimal.Decimal('902.346937180'),
           decimal.Decimal('920.012018681'), decimal.Decimal('937.332268954'),
           decimal.Decimal('954.498137713'), decimal.Decimal('971.926654816'),
           decimal.Decimal('989.451752425'), decimal.Decimal('1007.108936787'),
           decimal.Decimal('1037.157460451'), decimal.Decimal('1054.382933617'),
           decimal.Decimal('1073.052363158'), decimal.Decimal('1090.641976118'),
           decimal.Decimal('1108.038257122'), decimal.Decimal('1125.435575247'),
           decimal.Decimal('1143.022363425'), decimal.Decimal('1160.415113688'),
           decimal.Decimal('1178.040146590'), decimal.Decimal('1195.504595280'),
           decimal.Decimal('1213.258957863'), decimal.Decimal('1230.679465056')]




ft_diff = []
for i in range (0, len(ft_time) - 1):
  ft_diff.append(ft_time[i + 1] - ft_time[i])

print "Average duration: " + str(numpy.mean(ft_diff))
print "Iteration Number: " + str(len(ft_diff))

           
failure = 51
checkpoint = 44
shift = failure - checkpoint


after_ft_diff = []
for i in range (0, len(ft_time) - shift - 1):
  after_ft_diff.append(ft_time[i + shift + 1] - ft_time[i + shift])





Legends = []
Legends.append('Before Failure')
Legends.append('After Failure')

Parts = []

# line = plt.plot(ft_time[0:length], iter_num, color='g', linewidth=6)
# Parts.append(line[0])
# # line = plt.plot(normal_time[0:length], iter_num, '--', color='r', linewidth=1)
# # Parts.append(line[0])
# 
# plt.xlabel('Time (minute)', size=40)
# plt.ylabel('Iteration Number', size=40)
# plt.xticks(fontsize=30)
# plt.yticks(fontsize=30)
# # plt.axis([40, 160, 0, 0.03])
# plt.grid(True)
# plt.legend(Parts, Legends, loc='lower right', prop={'size':40})
# 
# plt.show()
# plt.close()

Parts = []

width = .5

line = plt.bar(range(0, failure), ft_diff[0:failure], width, color='k')
Parts.append(line[0])

new_iter_num = range(checkpoint, len(ft_diff)-shift)
for i in range(0, len(new_iter_num) - 1):
  new_iter_num[i] = new_iter_num[i] - width / 2.
line = plt.bar(new_iter_num, ft_diff[failure:len(ft_diff)], width, color='w')
Parts.append(line[0])



plt.xlabel('Main Loop Iteration Number')
plt.ylabel('Main Loop Duration (second)')
plt.xticks()
plt.yticks()
# plt.axis([40, 160, 0, 0.03])
plt.grid(True, axis='y')
plt.legend(Parts, Legends, loc='upper left', fontsize='small', frameon=False)

plt.annotate('frame 2', xy=(31, 20), xytext=(33, 40),
                arrowprops=dict(facecolor='black', shrink=0.05),
             fontsize='small')

plt.annotate('frame 3', xy=(54, 30), xytext=(48, 50),
                arrowprops=dict(facecolor='black', shrink=0.05),
             fontsize='small')

plt.annotate('checkpointing', xy=(44.5, 31), xytext=(46, 70),
                arrowprops=dict(facecolor='black', shrink=0.05),
             fontsize='small')

plt.annotate('rewinding', xy=(38, 100), fontsize='small')
# plt.annotate('rewinding', xy=(43.5, 85), xytext=(30, 115),
#                 arrowprops=dict(facecolor='black', shrink=0.05),
#              fontsize='small')

plt.annotate('failure', xy=(51, 13), xytext=(47, 30),
                arrowprops=dict(facecolor='black', shrink=0.05),
             fontsize='small')


plt.tight_layout(pad=0.1)



plt.show()


