#!/usr/bin/env python 

from optparse import OptionParser
import pprint
import sys

###############################################################################
#                               PARSER OPTIONS                                #
###############################################################################

parser = OptionParser()
parser.add_option('-i', '--in', dest='input', help='file containing data to parse')
parser.add_option('-o', '--out', dest='output', help='file to store results')

(options, args) = parser.parse_args()

###############################################################################
#                                    PARSE                                    #
###############################################################################

with open(options.input) as data:
 num_lines = len(data.readlines())

data = open(options.input)

def new_thread_table():
 table = {
   'started' : dict(),
   'delta_times' : dict(),
   'start_times' : dict()
   }
 return table

# dictionary to store all times and state
in_comp = dict()
copy_t = dict()
comp_t = dict()

# events to mark as beginning of copy jobs
comp_events = {'GAV stage', 'GAS stage', 'WIV stage', 'WIS stage'}

# parse each line, make state transitions and store times
num_line = 1
percent_check = 2
percent_inc = 2
sys.stdout.write("Parsed file % :   0")
sys.stdout.flush()
for line in data:
 if "wfcsize" in line or "rtcsize" in line:
  continue
 words = line.split(";")
 thread = words[0].strip()
 status = words[1].strip()
 event = words[2].strip()
 time = float(words[3].strip())
 if status == 'start':
  if event in comp_events:
   in_comp[thread] = True
  if thread in in_comp.keys() and in_comp[thread] == True:
   if thread not in comp_t.keys():
    comp_t[thread] = new_thread_table()
   table = comp_t[thread]
  else:
   if thread not in copy_t.keys():
    copy_t[thread] = new_thread_table()
   table = copy_t[thread]
   event = event + " " + " Copy job"
  table['started'][event] = True
  table['start_times'][event] = time
 else:
  if status == 'end':
   if thread in in_comp.keys() and in_comp[thread] == True:
    table = comp_t[thread]
   else:
    table = copy_t[thread]
    event = event + " " + " Copy job"
   assert(table['started'][event])
   if event not in table['delta_times'].keys():
    table['delta_times'][event] = 0
   table['delta_times'][event] = table['delta_times'][event] + (time - table['start_times'][event])
   if event in comp_events:
    in_comp[thread] = False
 num_line += 1
 if num_line * 100 / num_lines == percent_check:
  sys.stdout.write("\b\b\b%3i" % percent_check)
  sys.stdout.flush()
  percent_check += percent_inc
sys.stdout.write("\n")
sys.stdout.flush()

# save parsed times for different states
with open(options.output, 'w') as out:
 copy = {}
 comp = {}
 for thread in copy_t.keys():
  copy_pt = copy_t[thread]['delta_times']
  for event in copy_pt.keys():
   if event in copy.keys():
    copy[event] = copy[event] + copy_pt[event]
   else:
    copy[event] = copy_pt[event]
 for key in comp_t.keys():
  comp_pt = comp_t[key]['delta_times']
  for event in comp_pt.keys():
   if event in comp.keys():
    comp[event] = comp[event] + comp_pt[event]
   else:
    comp[event] = comp_pt[event]
 pprint.pprint('Copy Jobs:', stream=out)
 pprint.pprint(copy, width=1, stream=out)
 pprint.pprint('Compute Jobs:', stream=out)
 pprint.pprint(comp, width=1, stream=out)
