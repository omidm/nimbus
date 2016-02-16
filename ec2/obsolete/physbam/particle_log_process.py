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
particle_substep = [
	elem for elem in root.iter('scope')
	if elem.get('name') == 'Before stepping particle.']
print '#substep:', len(particle_substep)
time_array =[]
flag = False
temp = 0
for elem in particle_substep:
	if flag:
		time_array.append(temp+float(elem.find('time').get('value')))
		flag = False
		temp = 0
	else:
		temp = float(elem.find('time').get('value'))
		flag = True

print 'mean', numpy.mean(time_array[10:-10]), 'std', numpy.std(time_array[10:-10])
plt.plot(range(0, len(time_array)), time_array)
#plt.gca().set_ylim(ymin=0, ymax=0.08)
plt.show()
#plt.savefig('pic/figure'+sys.argv[1]+'.pdf')
