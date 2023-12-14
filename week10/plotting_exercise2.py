#!/usr/bin/env python

import numpy as np
import pandas as pd
from pydeseq2 import preprocessing
from matplotlib import pyplot as plt

art_df = pd.read_csv("artists.csv")
es_df = art_df[art_df['artist_nationality'] == 'Spanish']
es_list = es_df['artist_name'].value_counts().index

nats = art_df['artist_nationality'].value_counts().index
#print(nats)

lat_nat_list = ['Mexican', 'Cuban', 'Cuban-American', 'Uruguayan', 'Brazilian', 'Argentine', 'Columbian', 'Peruvian']

def count_lats(group):
    x = group['artist_nationality'][group['artist_nationality'].isin(lat_nat_list)].value_counts()
    return np.sum(x)

lat_count = art_df.groupby('year').apply(count_lats)
sizes = art_df.groupby('year')['artist_name'].apply(lambda x: len(x))
lat_share = lat_count / sizes

years = art_df['year'].value_counts().index
y_sort = sorted(years)

mx_df = art_df[art_df['artist_nationality'] == 'Mexican']

mx_count = mx_df.groupby('year')['artist_name'].apply(lambda x: len(x))
mx_share = mx_count / sizes

#lat_df = art_df['artist_nationality'][art_df['artist_nationality'].isin(lat_nat_list)]
#lat_df = art_df[art_df['artist_nationality'].isin(lat_nat_list)]
#lat_race = lat_df['artist_race'].to_numpy()

fig, ax = plt.subplots()
ax.plot(y_sort, lat_share, label = 'Latin Artists')
ax.plot(y_sort, mx_share, label = 'Mexican Artists')
ax.set_ylim(0, 0.25)
ax.set_ylabel('Proportion of Artist Entries')
ax.set_xlabel('Publication Year')
ax.legend()
ax.set_title('Proportion of Artists Featured in Janson and Gardner\n Art History Books who are of Latin Origin')
fig.savefig('Latin_artists.png')

fig2, ax2 = plt.subplots()
ax2.hist(sizes, bins = 6)
ax2.set_ylabel('Frequency')
ax2.set_xlabel('Number of Artist Entries')
ax2.set_title('Distribution of Entries in Janson and Gardner Art History Books')
fig2.savefig('entry_hist.png')



race_counts = art_df['artist_race'].value_counts()
races = art_df['artist_race'].value_counts().index
print(races)
xlabs = ['White', 'Black or\n African American', 'Asian', 'Native Hawaiian or\n Pacific Islander', 'American Indian or\n Alaska Native']
colors = []

fig3, ax3 = plt.subplots()
plt.bar(races, race_counts, color = 'black')
ax3.set_xticks(ticks = races, labels = xlabs)
ax3.set_ylabel('Number of Artist Entries')
ax3.set_xlabel('Artist Race')
ax3.set_title('Number of Artist Entries by Artist Race')
fig3.savefig('race_bar.png')

plt.show()
