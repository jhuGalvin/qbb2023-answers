#!/usr/bin/env python
import sys

fs = open(sys.argv[1])
lines = fs.readlines()

def findAvg(lines):
	nSum = 0

	for line in lines:
		num = int(line)
		nSum += num

	return nSum/len(lines)

print(findAvg(lines))