#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scipy.stats as sps
import statsmodels.formula.api as smf
import statsmodels.api as sm

# Exercise 1
# Wrangle the data, yeehaw
# Step 1.1
df = pd.read_csv("~/qbb2023-answers/week1/aau1043_dnm.csv")

# Step 1.2
sortDf = df.sort_values(by = ['Proband_id'])
dnms = df.groupby(['Proband_id', 'Phase_combined']).size().reset_index(name = 'Count')


deNovoCount = {}
for i in range(0, len(dnms), 2):
	deNovoCount[dnms.loc[i, 'Proband_id']] = [dnms.loc[i + 1, 'Count'], dnms.loc[i, 'Count']]

print(deNovoCount)

# Step 1.3
deNovoCountDF = pd.DataFrame.from_dict(deNovoCount, orient = 'index', columns = ['maternal_dnm', 'paternal_dnm'])

# Step 1.4
df_2 = pd.read_csv("~/qbb2023-answers/week1/aau1043_parental_age.csv", index_col = 'Proband_id')
print(df_2)

# Step 1.5
joinDf = pd.concat([deNovoCountDF, df_2], join = 'outer', axis = 1)
print(joinDf)

# Exercise 2
# Fit and interpret linear regression models with Python

# Step 2.1
fig, ax = plt.subplots()

#ax.scatter(joinDf["Mother_age"], joinDf['maternal_dnm'], label = "female")
#ax.set_ylabel("de Novo Mutations")
#ax.set_xlabel("Age")
#ax.set_title("de Novo Mutation Frequency by Maternal Age")
#fig.savefig("ex2_a.png")

ax.scatter(joinDf["Father_age"], joinDf['paternal_dnm'], label = "male")
ax.set_ylabel("de Novo Mutations")
ax.set_xlabel("Age")
ax.set_title("de Novo Mutation Frequency by Paternal Age")
fig.savefig("ex2_b.png")

# Step 2.2
matModel = smf.ols(formula = "maternal_dnm ~ 1 + Mother_age", data = joinDf)
matResults = matModel.fit()
print(matResults.summary())
print(matResults.pvalues)

# Step 2.3
pModel = smf.ols(formula = "paternal_dnm ~ 1 + Father_age", data = joinDf)
pResults = pModel.fit()
print(pResults.summary())
print(pResults.pvalues)

# Step 2.4 in readme

# Step 2.5
fig, ax = plt.subplots()
ax.hist(dnms.loc[dnms['Phase_combined'] == 'father', 'Count'], bins = 20, label = "Paternal", alpha = 0.5)
ax.hist(dnms.loc[dnms['Phase_combined'] == 'mother', 'Count'], bins = 20, label = "Maternal", alpha = 0.5)
ax.set_xlabel("de Novo Mutations")
ax.set_ylabel("Frequency")
ax.set_title('de Novo Mutations per Proband ID')
fig.savefig('ex2_c.png')

# Step 2.6
pData = dnms.loc[dnms['Phase_combined'] == 'father', 'Count']
mData = dnms.loc[dnms['Phase_combined'] == 'mother', 'Count']
print(sps.ttest_rel(pData, mData))
#print(sps.ttest_ind(pData, mData))




