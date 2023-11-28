#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 1.2 Plot top 2 genotype PCs
PC1 = []
PC2 = []

for line in open('genotypePCs.eigenvec'):
	line = line.rstrip('\n').split()

	PC1.append(float(line[2]))
	PC2.append(float(line[3]))

fig, ax = plt.subplots()
ax.scatter(PC1, PC2)
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("Top 2 PCs for genotypes.vcf")

plt.show()
fig.savefig("genotypesVcf_PCs")

MAFs = []

# 2.2 Plot AFs
f = open('plink.frq', 'r')

for line in f:
	line = line.split()
	if 'MAF' in line:
		continue

	MAFs.append(line[4])

MAFs = np.array(MAFs).astype(float)

counts, bins = np.histogram(MAFs, bins = 20)
afHist = plt.stairs(counts, bins)
plt.title('Sample Allele Frequency Spectrum')
plt.xlabel('Allele Frequency')
plt.savefig('AFs')


# Function for Step 3.2
def makeManhat(file):
	gwasResults = []
	headers = []

	with open(file, 'r') as f:
		for line in f:
			line = line.split()
			if 'CHR' in line:
				headers = [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8],]

			if 'CHR' not in line:
				row = [int(line[0]), line[1], int(line[2]), line[3], line[4], line[5], float(line[6]), float(line[7]), float(line[8])]
				gwasResults.append(row)


	del gwasResults[0]
	gwasDf = pd.DataFrame(gwasResults, columns = headers)
	del gwasResults

	# Create column of -log10 p values
	gwasDf['-log10P'] = -np.log10(gwasDf['P'])
	gwasDf.CHR = gwasDf.CHR.astype('category')
	gwasDf = gwasDf.sort_values('CHR', ascending = True)

	#
	gwasDf['ind'] = range(len(gwasDf))
	gwas_Grouped = gwasDf.groupby('CHR')

	filtDf = gwasDf[gwasDf['-log10P'] >= 4.0]
	filtDf_Grouped = filtDf.groupby('CHR')

	# Make Manhattan plot
	fig = plt.figure(figsize = (14, 8)) # Set the figure size
	ax = fig.add_subplot(111)
	colors = ['lightsteelblue', 'steelblue']
	x_labels = []
	x_labels_pos = []

	for num, (name, group) in enumerate(gwas_Grouped):
		group.plot(kind = 'scatter', x = 'ind', y = '-log10P', color = colors[num % len(colors)], ax = ax)
		x_labels.append(name)
		x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))

	for num, (name, group) in enumerate(filtDf_Grouped):
		group.plot(kind = 'scatter', x = 'ind', y = '-log10P', color = 'orange', ax = ax)

	ax.set_xticks(x_labels_pos)
	ax.set_xticklabels(x_labels)

	# Axis limits
	ax.set_xlim([0, len(gwasDf)])
	ax.set_ylim([0, 10])

	# Axis label
	ax.set_xlabel('Chromosome')

	# Set title
	title = file.split('_gwas_results.assoc.linear')[0]
	ax.set_title(title)
	fig.savefig(title)


# Function for Step 3.3
def findSNP(assocFile, vcfFile):
	gwasResults = []
	headers = []

	with open(assocFile, 'r') as f:
		for line in f:
			line = line.split()
			if 'CHR' in line:
				headers = [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8],]

			if 'CHR' not in line:
				row = [int(line[0]), line[1], int(line[2]), line[3], line[4], line[5], float(line[6]), float(line[7]), float(line[8])]
				gwasResults.append(row)

	del gwasResults[0]
	gwasDf = pd.DataFrame(gwasResults, columns = headers)
	del gwasResults

	gwasDf['-log10P'] = -np.log10(gwasDf['P'])
	gwasDf = gwasDf.sort_values('-log10P', ascending = True)
	maxP = gwasDf.iloc[len(gwasDf) - 1]
	maxSNP = maxP.SNP
	maxSNPs = gwasDf[gwasDf['SNP'].str.contains(maxSNP)]
	#print(maxSNPs)

	GTList = []
	snpList = []

	for line in open(vcfFile, 'r'):
		if line.startswith('#'):
			continue
		if maxSNP in line:
			line = line.rstrip('\n').split('\t')
			GTList = line

	for i in range(9, len(GTList)):
		if GTList[i] == './.':
			snpList.append('NA')
		if GTList[i] == '0/0':
			snpList.append('AA')
		if GTList[i] == '0/1':
			snpList.append('AG')
		if GTList[i] == '1/1':
			snpList.append('GG')
	#print(snpList)

	pheno = pd.read_csv('CB1908_IC50.txt', sep = '\t')
	pheno['GT'] = snpList

	AAdf = pheno[pheno['GT'] == 'AA']
	AGdf = pheno[pheno['GT'] == 'AG']
	GGdf = pheno[pheno['GT'] == 'GG']

	AAval = AAdf['CB1908_IC50'].values
	AGval = AGdf['CB1908_IC50'].values
	GGval = GGdf['CB1908_IC50'].values

	nAGval = np.delete(AGval, [12, 30])
	#print(nAGval)

	labels = ['AA', 'AG', 'GG']
	values = [AAval, nAGval, GGval]
	plt.boxplot(values, vert = True, labels = labels)
	plt.title('IC 50 Scores by rs10876043 Genotype')
	plt.xlabel('Genotype')
	plt.ylabel('IC 50')

	plt.savefig('rs10876043')

# Step 3.2 Visualizing GWAS results
makeManhat('GS451_IC50_gwas_results.assoc.linear')
makeManhat('CB1908_IC50_gwas_results.assoc.linear')

# Couldnt get them to load as panels in the same figure, but hopefully they still look pretty enough!
#fig, (ax1, ax2) = plt.subplots(2)
#fig.suptitle('GWAS Results')
#plt.plot(ax1, ax2)
#plt.show()

# Step 3.3 Visualizing effect-size
findSNP('CB1908_IC50_gwas_results.assoc.linear', 'genotypes.vcf')
