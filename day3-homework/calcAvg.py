#!/usr/bin/env python


testSet = [1, 2, 6, 4, 5]

def findAvg(numSet):
	nSum = 0
	for i in numSet:
		nSum += i

	return nSum/len(numSet)

print(findAvg(testSet))