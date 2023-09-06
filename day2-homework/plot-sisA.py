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

# Iterate through
expr = []
fExpr = []
mExpr = []

for i in range(len(samples)):
    expr.append(data[row,i])
    if "female" in samples[i]:
        fExpr.append(expr[i])
    else:
        mExpr.append(expr[i])

devStages = ['10', '11', '12', '13', '14A', '14B', '14C', '14D']

# Plot data
fig, ax = plt.subplots()
ax.set_title( "SisA (FBtr0073461)" )
ax.set_ylabel("mRNA abundance (RPKM)")
ax.set_xlabel("developmental stage")

ax.plot(devStages, fExpr)
ax.plot(devStages, mExpr)

plt.xticks(rotation = 90)
plt.tight_layout()
fig.savefig( "FBtr0073461.png" )
plt.close( fig )