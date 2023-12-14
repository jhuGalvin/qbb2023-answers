#!/usr/bin/env python

import numpy as np
import pandas as pd
from pydeseq2 import preprocessing
from matplotlib import pyplot as plt

# read in data
counts_df = pd.read_csv("gtex_whole_blood_counts_formatted.txt", index_col = 0)

# read in metadata
metadata = pd.read_csv("gtex_metadata.txt", index_col = 0)

# normalize
counts_df_normed = preprocessing.deseq2_norm(counts_df)[0]

# log
counts_df_logged = np.log2(counts_df_normed + 1)

# merge with metadata
full_design_df = pd.concat([counts_df_logged, metadata], axis = 1)

# Step 1.1
subjDf = full_design_df.loc['GTEX-113JC']
subjDf = subjDf[subjDf != 0]
subjDf = subjDf.iloc[:-3]

fig1, ax1 = plt.subplots()
ax1.hist(subjDf, bins = 20)
ax1.set_xlabel('Gene Expression, Normalized')
ax1.set_ylabel('Frequency')
ax1.set_title('Total Gene Expression in GTEX Subject 113JC')
fig1.savefig('GTEX_113JC.png')

# Step 1.2
maleDf = full_design_df.loc[full_design_df['SEX'] == 1, 'MXD4']
femaleDf = full_design_df.loc[full_design_df['SEX'] == 2, 'MXD4']

fig2, ax2 = plt.subplots()
ax2.hist(maleDf, bins = 15, label = 'Male')
ax2.hist(femaleDf, bins = 15, label = 'Female')
ax2.set_xlabel('Gene Expression, Normalized')
ax2.set_ylabel('Frequency')
ax2.set_title('MXD4 Expression Stratified by Sex')
ax2.legend()
fig2.savefig('MDX4.png')

# Step 1.3
ages = full_design_df['AGE'].value_counts()
age_counts = ages.values
ages = ages.index
sorted_ages = [ages[3], ages[4], ages[2], ages[1], ages[0], ages[5]]

fig3, ax3 = plt.subplots()
ax3.bar(sorted_ages, height = age_counts)
ax3.set_xlabel('Age Ranges')
ax3.set_ylabel('Frequency')
ax3.set_title('Number of Subjects per Age Group')
fig3.savefig('Age_Dist.png')

# Step 1.4
LPXN_male_df = full_design_df[full_design_df['SEX'] == 1]
m_med = LPXN_male_df.groupby('AGE')['LPXN'].apply(lambda x: x.median())

LPXN_female_df = full_design_df[full_design_df['SEX'] == 2]
f_med = LPXN_female_df.groupby('AGE')['LPXN'].apply(lambda x: x.median())

fig4, ax4 = plt.subplots()
ax4.plot(sorted_ages, m_med, label = 'Male')
ax4.plot(sorted_ages, f_med, label = 'Female')
ax4.set_xlabel('Age Groups')
ax4.set_ylabel('Median LPXN Expression')
ax4.set_title('Median LPXN Expression across Age Groups')
ax4.legend()

fig4.savefig('LPXN_med.png')

# Exercise 2: Independent data visualization





