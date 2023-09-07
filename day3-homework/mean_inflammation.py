#!/usr/bin/env python

import sys

p = sys.argv[1]
f = open(sys.argv[2])
lines = f.readlines()

def meanInflammation(patient, file):
	patient = int(patient)
	patientData = lines[patient].split(",")
	iSum = 0
	for i in range(len(patientData)):
		iSum += float(patientData[i])
	iAvg = iSum/len(patientData)
	return iAvg

print(meanInflammation(p,f))