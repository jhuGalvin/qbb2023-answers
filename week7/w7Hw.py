#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

# Part 3a
ONT = []
bisulph = []

for line in open("ONT.cpg.chr2.bedgraph"):
	data = line.rstrip().split()
	ONT.append(float(data[1]))

for line in open("bisulfite.cpg.chr2.bedgraph"):
	data = line.rstrip().split()
	bisulph.append(float(data[1]))

ONT_set = set()
for i in range(len(ONT)):
    if ONT[i] not in ONT_set:
        ONT_set.add(ONT[i])

bisulph_set = set()
for i in range(len(bisulph)):
    if bisulph[i] not in bisulph_set:
        bisulph_set.add(bisulph[i])

shared = ONT_set.intersection(bisulph_set)

# Part 3b
ONT_coverage = []
bisulph_coverage = []

for line in open("ONT.cpg.chr2.bedgraph"):
	data = line.rstrip().split()
	ONT_coverage.append(float(data[4]))

for line in open("bisulfite.cpg.chr2.bedgraph"):
	data = line.rstrip().split()
	bisulph_coverage.append(float(data[4]))

fig, ax1 = plt.subplots()
ax1.hist(ONT_coverage, bins = 2000, color = 'blue', alpha = 0.5, label = "ONT Coverage")
ax1.set_xlabel("Read Coverage")
ax1.set_ylabel("Frequency")
ax1.hist(bisulph_coverage, bins = 2000, color = 'orange', alpha = 0.5, label = "Bisulphite Coverage")
ax1.set_title("Coverage Distribution of ONT and Bisulphite Reads")
ax1.set_xlim(0,100)
ax1.legend()

fig.savefig("method_coverage.png")
fig.tight_layout()

# Part 3c
ONTme = []
bisulphMe = []
for line in open("ONT.cpg.chr2.bedgraph"):
	data = line.rstrip().split()
	ONTme.append(float(data[3]))

for line in open("bisulfite.cpg.chr2.bedgraph"):
	data = line.rstrip().split()
	bisulphMe.append(float(data[3]))

ONT_score = []
bisulph_score = []

for i in range(len(ONT)):
	if ONT[i] in shared:
		ONT_score.append(ONTme[i])

for j in range(len(bisulph)):
	if bisulph[j] in shared:
		bisulph_score.append(bisulphMe[j])

mScoreCorr = np.corrcoef(ONT_score, bisulph_score)

fig = plt.figure()
plt.imshow(np.log10(1 + np.histogram2d(ONT_score, bisulph_score, bins = 100)[0]))
plt.xlabel("ONT Methylation Score")
plt.ylabel("Bisulphsulphite Methylation Score")
plt.title(f'Correlation of ONT and Bisulphite Methylation Scores, R = {mScoreCorr[0,1]:0.3f}')
plt.show()
fig.savefig("ONT_bisulphite_Corr.png")


# Part 3d
ONT_norm = []
ONT_tumor = []
ONT_cov = []
ONT_tumorCov = []
ONT_score = []
ONT_tumorScore = []

for line in open("normal.ONT.chr2.bedgraph"):
	data = line.rstrip().split()
	ONT_norm.append(float(data[1]))
	ONT_score.append(float(data[3]))
	ONT_cov.append(int(data[4]))

for line in open("tumor.ONT.chr2.bedgraph"):
	data = line.rstrip().split()
	ONT_tumor.append(float(data[1]))
	ONT_tumorScore.append(float(data[3]))
	ONT_tumorCov.append(int(data[4]))

ONT_set = set()

for i in range(len(ONT_norm)):
    if ONT_norm[i] not in ONT_set:
        ONT_set.add(ONT_norm[i])


ONT_tumorSet = set()

for i in range(len(ONT_tumor)):
    if ONT_tumor[i] not in ONT_tumorSet:
        ONT_tumorSet.add(ONT_tumor[i])

bisulph_shared = ONT_set.intersection(ONT_tumorSet)

ONT_mNorm = []
ONT_mTumor = []
for i in range(len(ONT_norm)):
	if ONT_norm[i] in bisulph_shared:
		ONT_mNorm.append(ONT_score[i])

for j in range(len(ONT_tumor)):
	if ONT_tumor[j] in bisulph_shared:
		ONT_mTumor.append(ONT_tumorScore[j])

total_ONT_mChange = []

for i in range(len(ONT_mNorm)):
	ONT_mChange = ONT_mNorm[i] - ONT_mTumor[i]
	if ONT_mChange != 0:
		total_ONT_mChange.append(ONT_mChange)

bisulph_norm = []
bisulph_tumor = []
bisulph_cov = []
bisulph_tumorCov = []
bisulph_score = []
bisulph_tumorScore = []


for line in open("normal.bisulfite.chr2.bedgraph"):
	data = line.rstrip().split()
	bisulph_norm.append(float(data[1]))
	bisulph_score.append(float(data[3]))
	bisulph_cov.append(int(data[4]))

for line in open("tumor.bisulfite.chr2.bedgraph"):
	data = line.rstrip().split()
	bisulph_tumor.append(float(data[1]))
	bisulph_tumorScore.append(float(data[3]))
	bisulph_tumorCov.append(int(data[4]))

bisulph_set = set()
for i in range(len(bisulph_norm)):
    if bisulph_norm[i] not in bisulph_set:
        bisulph_set.add(bisulph_norm[i])


bisulph_tumorSet = set()

for i in range(len(bisulph_tumor)):
    if bisulph_tumor[i] not in bisulph_tumorSet:
        bisulph_tumorSet.add(bisulph_tumor[i])

bisulph_shared = bisulph_set.intersection(bisulph_tumorSet)

bisulph_mNorm = []
bisulph_mTumor = []

for i in range(len(bisulph_norm)):
	if bisulph_norm[i] in bisulph_shared:
		bisulph_mNorm.append(bisulph_score[i])

for j in range(len(bisulph_tumor)):
	if bisulph_tumor[j] in bisulph_shared:
		bisulph_mTumor.append(bisulph_tumorScore[j])

total_bisulph_mChange = []
for i in range(len(bisulph_mNorm)):
	bisulph_mChange = bisulph_mNorm[i] - bisulph_mTumor[i]
	if bisulph_mChange != 0:
		total_bisulph_mChange.append(bisulph_mChange)

mChangeCorr = np.corrcoef(total_ONT_mChange, total_bisulph_mChange[0:len(total_ONT_mChange)])

fig, ax1 = plt.subplots()
ax1.violinplot([total_ONT_mChange, total_bisulph_mChange])
ax1.set_xticks([1, 2])
ax1.set_xticklabels(['ONT', 'Bisulphite'])
ax1.set_xlabel("Sequencing Method")
ax1.set_ylabel("Methylation Change (Tumor - Normal)")
ax1.set_title(f'Distribution of Methylation Changes, R = {mChangeCorr[0,1]:0.3f}')
fig.savefig("mChange_violin.png") 
fig.tight_layout()
plt.close
