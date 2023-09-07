#!/usr/bin/env python

import sys

p1 = sys.argv[1]
p2 = sys.argv[2]
f = open(sys.argv[3])
lines = f.readlines()

def meanInflammation(patient, file):
	patient = int(patient)
	patientData = lines[patient].split(",")
	iSum = 0
	for i in range(len(patientData)):
		iSum += float(patientData[i])
	iAvg = iSum/len(patientData)
	print("Patient", patient, "average inflammation: ", iAvg)

def differenceInflammation(patient1, patient2, file):
	patient1 = int(patient1)
	patient2 = int(patient2)

	patient1Data = lines[patient1].split(",")
	patient2Data = lines[patient2].split(",")

	for i in range(len(patient1Data)):
		diff = abs(int(patient1Data[i]) - int(patient2Data[i]))
		print("Day", i + 1, ": ", diff)


meanInflammation(p1,f)
differenceInflammation(p1, p2, f)