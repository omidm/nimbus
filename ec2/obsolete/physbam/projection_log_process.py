#!/usr/bin/env python

import xml.etree.ElementTree as ET
import numpy
import matplotlib.pyplot as plt
import sys

def ExtractStep(elem, name):
	return numpy.array([float(step.find('time').get('value'))
		for step in elem.iter('scope')
		if step.get('name') == name])
	

tree = ET.parse('temp/log_'+sys.argv[1]+'.txt')
root = tree.getroot()
projection_substep = [
	elem for elem in root.iter('scope')
	if elem.get('name') == 'Implicit Part']
print '#substep:', len(projection_substep)
mean_array =[]
std_array =[]
for elem in projection_substep:
	execute_1 = ExtractStep(elem, 'Enter iteration.')
	execute_2 = ExtractStep(elem, 'After second barrier.')
	execute_3 = ExtractStep(elem, 'After third barrier.')
	execute_4 = ExtractStep(elem, 'After fourth barrier.')
	execute_5 = ExtractStep(elem, 'After sixth barrier.')
	execute_all = execute_1+execute_2+execute_3+execute_4+execute_5
	mean_array.extend(execute_all[3:-3])
	#mean_array.append(numpy.mean(execute_all[3:-3]))
	#std_array.append(numpy.std(execute_all[3:-3]))

#print mean_array
#print std_array
print 'mean', numpy.mean(mean_array), 'std', numpy.std(mean_array)
#plt.errorbar(range(0, len(mean_array)), mean_array, std_array)
plt.plot(range(0, len(mean_array)), mean_array)
plt.gca().set_ylim(ymin=0, ymax=0.08)
plt.show()
#plt.savefig('pic/figure'+sys.argv[1]+'.pdf')
