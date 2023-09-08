#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Get starting frequency and population size
# Input parameters for function

# make list to store allele frequencies

# Return list of allele frequencies at each timepoint
# Stop when allele reaches fixation, return timepoint (length of list)

# while allele frequency is between 0 & 1:
# 	Get new allele frequency for next generation
# 	by drawing from binomial distribution;
# 	(convert number of successes into a frequency);
#
# 	store allele frequencies in allele frequency list

def drift(af, pop):

	afGen = []

	while 0 < af < 1:
		success = np.random.binomial(2*pop, af)
		af = success/(2*pop)
		afGen.append(af)

	nGen = len(afGen)
	return afGen, nGen

results = drift(0.3, 300)
print(results)

fig, ax = plt.subplots()
for i in range(40):
	results = drift(0.3, 300)
	x = range(0, results[-1])
	y = results[0]
	ax.plot(x, y)

fig.savefig("afModel.png")
ax.set_ylabel("Allele Frequency")
ax.set_xlabel("Generations")
ax.set_title("Time to Fixation of Allele")
plt.show()




#plt.show()