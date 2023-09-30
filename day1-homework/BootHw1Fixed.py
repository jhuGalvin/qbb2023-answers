#!/usr/bin/env python

import numpy

# Exercise 1
f = open("inflammation-01.csv", "r")
lines = f.readlines()
patient5 = lines[4].split(",")
print("First day:" + patient5[0] + " Tenth day:" + patient5[9] + " Last day:" + patient5[-1])

Ex2Avgs = []

# Exercise 2
print("Average flare-ups per patient:")
for i in range(10):
	patient = lines[i].split(",")
	for j in range(len(patient)):
		patient[j] = int(patient[j])
	patientAvg = numpy.average(patient)
	Ex2Avgs.append(patientAvg)
	print(patientAvg)

# Exercise 3
patientMax = max(Ex2Avgs)
patientMin = min(Ex2Avgs)
print("\nMinimum and maximum flare-ups: \n" + str(patientMin) + " " + str(patientMax))

# Exercise 4
patient1 = lines[0].split(",")
print("\nDifference in daily flare-ups between Patient 1 & Patient 5:")
#print(len(patient1) == len(patient5))
for i in range(len(patient1)):
	patient1[i] = int(patient1[i])
	patient5[i] = int(patient5[i])
	patientDiff = abs(patient1[i] - patient5[i])
	print(patientDiff)

f.close()