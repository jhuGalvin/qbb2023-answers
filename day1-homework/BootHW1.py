#!/usr/bin/env python

f = open("inflammation-01.csv", "r")
lines = f.readLines()

for line in lines:
	print(line)

f.close()