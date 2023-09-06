#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Get dataset to recreate Fig 3B from Lott et al 2011 PLoS Biology https://pubmed.gov/21346796
# wget https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/bulk_RNA-seq/extra_data/all_annotated.csv

transcripts = np.loadtxt( "all_annotated.csv", delimiter=",", usecols=0, dtype="<U30", skiprows=1 )
print( "transcripts: ", transcripts[0:5] )

samples = np.loadtxt( "all_annotated.csv", delimiter=",", max_rows=1, dtype="<U30" )[2:]
print( "samples: ", samples[0:5] )

data = np.loadtxt( "all_annotated.csv", delimiter=",", dtype=np.float32, skiprows=1, usecols=range(2, len(samples) + 2) )
print( "data: ", data[0:5, 0:5] )

# Find row with transcript of interest
for i in range(len(transcripts)):
    if transcripts[i] == 'FBtr0073461':
        row = i

# Find columns with female samples
fCols = []
for i in range(len(samples)):
    if "female" in samples[i]:
        fCols.append(i)

# Subset female data
fExpression = data[row, fCols]

# Prepare data
fX = samples[fCols]
fY = fExpression

# Find columns with male samples
mCols = []
for i in range(len(samples)):
    if "female" in samples[i]:
        mCols.append(i)

# Subset male data
mExpression = data[row, mCols]

# Prepare data
mX = samples[mCols]
mY = mExpression

# Plot data
fig, ax = plt.subplots()
ax.set_title( "SisA (FBtr0073461)" )
ax.set_ylabel("mRNA abundance (RPKM)")
ax.set_xlabel("developmental stage")

ax.plot(fX,fY)
ax.plot(mX,mY)

plt.xticks(rotation = 90)
#ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, ha = 'right')
plt.tight_layout()
fig.savefig( "FBtr0073461.png" )
plt.close( fig )