#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np

bait_file, washU_file, usc_file = sys.argv[1:]

totalData = []
baitData = []
baitCoords = []
baitNames = []

with open(bait_file, 'r') as bf:
	for line in bf:
		data = line.split()
		baitCoord = data[1] + ':' + data[2]
		baitName = data[4]
		baitCoords.append(baitCoord)
		baitNames.append(baitName)

bait_df = pd.DataFrame({'coords':baitCoords, 'names':baitNames})

with open(washU_file, 'r') as wf:
	for line in wf:
		data = line.split()
		value = data[2]
		inter1 = data[0].split(',')
		inter2 = data[1].split(',')
		#set1 = min(set1, set2)
		#set2 = max(set1, set2)

		# Produce the .bed file fields No. 1. to 8.
		interData = [inter1[0], min(int(inter1[1]), int(inter2[1])), max(int(inter1[2]), int(inter2[2])), '.', 'score', value, '.', 0]
		coords = [inter1[1], inter1[2]]
		coords2 = [inter2[1], inter2[2]]

		c_coord = inter1[1] + ':' + inter1[2]
		c_coord2 = inter2[1] + ':' + inter2[2]

		bf_data = []

		# Produce the .bed file fields No. 9. to 18.
		# (sourceChrom, sourceStart, sourceEnd, sourceName, sourceStrand, 
		# targetChrom, targetStart, targetEnd, targetName, targetStrand)

		if c_coord in baitCoords and c_coord2 in baitCoords:
			condition = bait_df['coords'] == c_coord
			name = bait_df.loc[condition]['names']
			name = name.to_string(header = False, index = False)

			condition2 = bait_df['coords'] == c_coord2
			name2 = bait_df.loc[condition2]['names']
			name2 = name2.to_string(header = False, index = False)

			bf_data = [interData[0], coords[0], coords[1], name, '+', interData[0], coords2[0], coords2[1], name2, '+']

		if c_coord not in baitCoords and c_coord2 in baitCoords:
			condition = bait_df['coords'] == c_coord2
			name = bait_df.loc[condition]['names']
			name = name.to_string(header = False, index = False)

			bf_data = [interData[0], coords[0], coords[1], name, '+', interData[0], coords2[0], coords2[1], '.', '-']

		if c_coord in baitCoords and c_coord2 not in baitCoords:
			condition = bait_df['coords'] == c_coord
			name = bait_df.loc[condition]['names']
			name = name.to_string(header = False, index = False)

			bf_data = [interData[0], coords[0], coords[1], name, '+', interData[0], coords2[0], coords2[1], '.', '-']

		combined_data = interData + bf_data
		
		totalData.append(combined_data)

cols = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'value', 'exp', 'color', 'sourceChrom', 'sourceStart', 
		 'sourceEnd', 'sourceName', 'sourceStrand', 'targetChrom', 'targetStart', 'targetEnd', 'targetName', 'targetStrand']

df = pd.DataFrame(totalData, columns = cols)

df['value'] = df['value'].astype(float)
#df['max']= df.groupby('sourceName')['value'].transform('max')
maxVal = max(df['value'].value_counts())
print(maxVal)
df['score'] = (df['value'] / maxVal) * 1000
df['score'] = df['score'].astype(int)
df = df.copy()

#df.drop('max', axis = 1, inplace = True)

df.to_csv(usc_file, sep = '\t', header = False, index = False)

# Add header
header = 'track type=interact name="pCHIC" description="Chromatin interactions" useScore=on maxHeightPixels=200:100:50 visibility=full\n'

with open(usc_file, 'r') as f:
    ucsc_f = f.read()

header_ucsc = header + ucsc_f

# Write the updated content back to the file
with open(usc_file, 'w') as f:
    f.write(header_ucsc)

# Step 2.2
df2 = pd.read_csv('chicago_Results_UCSC.bed', names = cols, sep = '\t')

top_bb = df2[df2['targetStrand'] == '+']
top_bb = top_bb.copy()
top_bb.sort_values(by = 'score', ascending = False, inplace = True)

top_bt = df2[df2['targetStrand'] == '-']
top_bt = top_bt.copy()
top_bt.sort_values(by = 'score', ascending = False, inplace = True)

top_bb.to_csv('max_bb.txt', sep = ' ', header = True, index = False)
top_bt.to_csv('max_bt.txt', sep = ' ', header = True, index = False)




