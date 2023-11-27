#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Initialize empty lists
DPs = []
GQs = []
AFs = []
EFs = []

# Step 3.0 Parse vcf file
for line in open('annotVar.vcf', 'r'):

    if line.startswith('#'):
        continue

    fields = line.rstrip('\n').split('\t')

    stats = fields[7].split(';')
    gtStats = fields[14].split(':')
    EFHeader = stats[-1].split('|')


    # Grab what you need from `fields`
    # 3.1 Read depth (DP)
    DP = stats[7].split('DP=')[1]
    DPs.append(DP)

    # 3.2 Genotype quality (GQ)
    GQs.append(gtStats[1])

    # 3.3 Allele frequency (AF)
    AF = stats[3].split('AF=')[1]
    comma = ','

    if comma in AF:
    	continue
    else:
    	AFs.append(AF)

    # 3.4 Predicted Effect (EF) from snpEff
    EF = EFHeader[1]

    if '&' in EF:
    	dubEF = EF.split('&')
    	EFs.append(dubEF[0])
    	EFs.append(dubEF[1])
    else:
    	EFs.append(EF)

uGenVar = EFs.count('upstream_gene_variant')
dGenVar = EFs.count('downstream_gene_variant')
synVar = EFs.count('synonymous_variant')
misVar = EFs.count('missense_variant')
spliceVar = EFs.count('splice_region_variant')
sR = EFs.count('stop_retained_variant')
sLost = EFs.count('stop_lost')
sGained = EFs.count('stop_gained')
DisVar = EFs.count('disruptive_inframe_insertion')

snpEff_for_bar = [synVar, uGenVar, misVar, dGenVar, spliceVar, sR, sLost, sGained, DisVar]
snpEffNames = ['Synonymous\n Variant', 'Upstream\n Gene\n Variant', 'Missense\n Variant', 'Downstream\n Gene\n Variant', 'Splice\n Region\n Variant', 'Stop\n Retained', 'Stop\n Lost', 'Stop\n Gained', 'Disruptive\n Inframe\n Insertion']

#3.4 Graph summary of features & predicted effects
fig, ((DPax, GQax), (AFax, Eax)) = plt.subplots(2, 2)
fig.suptitle('Summary of Variant Features & Predicted Effects')

DPax.hist(DPs, bins = 20)
DPax.set_title('Variant Read Depth')
DPax.tick_params(axis = 'x', which = 'both', bottom = False, labelbottom = False)

GQax.hist(GQs, bins = 12)
GQax.set_title('Variant Genotype Quality')
GQax.tick_params(axis = 'x', which = 'both', bottom = False, labelbottom = False)

AFax.hist(AFs, bins = 12)
AFax.set_title('Variant Allele Frequency')
AFax.tick_params(axis = 'x', which = 'both', bottom = False, labelbottom = False)

Eax.bar(snpEffNames, snpEff_for_bar)
eBar = Eax.bar(snpEffNames, snpEff_for_bar)
Eax.bar_label(eBar, labels = snpEff_for_bar)
Eax.tick_params(axis = 'x', labelsize = 6)
Eax.set_title('Predicted Variant Effects')

plt.show()
fig.savefig("VariantGraph")



