# Assignment 7: Nanopore Sequencing and Methylation
## Part 2: Exploring CpG methylation
(Q1): A majority of the CpG islands are methylated, as the two tracks (CpG islands & bisulphite CpG) are nearly identical. 

## Part 3: Comparing nanopore vs. bisulphite sequencing methylation calling

Part 3b 
(Q2):The bisulphite sequencing has more reads with greater coverage, but ONT sequencing has more homogeneity in coverage (smaller spread). To me, the bisulphite sequencing looks like it is better to use, but may require some quality control before further analyses. Still, the ONT sequencing still appears good to use, especially considering ONT's cost, speed, and ease of use.

Part 3d 
(Q3): Peaks resemble one another (same locations & relative sizes in each datasets) between the two sequencing methods, but there is a greater number of detected methylation events in the ONT dataset. this may be due to the increased size of the ONT dataset relative to the bisulphite sequencing dataset(?). 

(Q4): There is both increased and decreased methylation in tumor samples, with two notable peaks on either side of the 0 methylation change value. Cancer is associated with differential methylation, but usually not as hypo- or hypermethylation exclusively. To understand the effects of differential methylation on tumorigenesis, we would need to look at the differentially methylated sites and determine their role (i.e. are they promoters?, enhancers?, genes?, etc.).

Part 4a 
(Q5): A region in Exon 1 of the DNMT3A gene is hypermethylated in cancer samples (visible in the snapshot, more blue in the normal sample, more red in the tumor sample), repressing the gene. Because this gene is a methyltransferase, silencing of this gene likely prevents CpG islands methylation in genes whose functions must be controlled throughout development. Unregulated expression of these genes is probably what is responsible for many cancers.  

Part 4b 
(Q6): An imprinted gene is expressed or silenced depending on which parent the allele came from, either partially or completely. Imprinted genes often have developmental roles. 

(Q7): When reads are phased, the reads are clustered in haplotypes based on variants which segregate together in reads. 

(Q8): Reads can only be clustered if they have enough variants, which means there must be enough structured variation to separate into haplotypes, at least two as we did in this assignment. At certain regions within Exon 1, there was not enough information for IGV to cluster the reads, mainly because there were not enough, or any, variants present.


