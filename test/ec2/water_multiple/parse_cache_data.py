#!/usr/bin/env python 

from optparse import OptionParser
import pprint

parser = OptionParser()
parser.add_option('-i', '--in', dest='input', help='file containing data to parse')
parser.add_option('-o', '--out', dest='output', help='file to store results')

(options, args) = parser.parse_args()

data = open(options.input)

def new_thread_table():
 table = {
   'started' : dict(),
   'delta_times' : dict(),
   'start_times' : dict()
   }
 return table

in_copy = dict()
copy_t = dict()
comp_t = dict()

copy_events = {'LC', 'RCS', 'RCR'}

active_copy = {}
for line in data:
 if "region" in line or "size" in line:
  continue
 words = line.split(";")
 thread = words[0].strip()
 status = words[1].strip()
 event = words[2].strip()
 time = float(words[3].strip())
 if status == 'start':
  if event in copy_events:
   in_copy[thread] = True
   active_copy[thread] = event
  if in_copy[thread] == True:
   if thread not in copy_t.keys():
    copy_t[thread] = new_thread_table()
   table = copy_t[thread]
   if event not in copy_events and active_copy[thread] and active_copy[thread] != None:
       event = event + " " + active_copy[thread]
  else:
   if thread not in comp_t.keys():
    comp_t[thread] = new_thread_table()
   table = comp_t[thread]
  table['started'][event] = True
  table['start_times'][event] = time
 else:
  if status == 'end':
   if in_copy[thread] == True:
    table = copy_t[thread]
    if event not in copy_events and active_copy[thread] and active_copy[thread] != None:
        event = event + " " + active_copy[thread]
   else:
    table = comp_t[thread]
   assert(table['started'][event])
   if event not in table['delta_times'].keys():
    table['delta_times'][event] = 0
   table['delta_times'][event] = table['delta_times'][event] + (time - table['start_times'][event])
   if event in copy_events:
    in_copy[thread] = False
    active_copy[thread] = None

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
