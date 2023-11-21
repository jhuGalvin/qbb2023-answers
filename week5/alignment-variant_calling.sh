#!/bin/bash

# Step 1.1 index reference genome
bwa index sacCer3.fa

# Steps 1.2 & 1.3, align reads to reference, sort and index sam/bam files
for sample in *.fastq
do 
	echo 'Aligning sample:' ${sample}
	bwa mem -t 4 -R "@RG\tID:${sample}\tSM:${sample}" \
		sacCer3.fa \
		${sample} -o ${sample}.sam
		samtools sort ${sample}.sam -o ${sample}.bam
		samtools index ${sample}.bam
done

# Exercise 2
# Step 2.1 Call variants
ls *.bam > yeastBam_list
freebayes -f sacCer3.fa -p 4 -= -L yeastBam_list > yeastVar.vcf

# Step 2.2 Filter vcf file
vcffilter yeastVar.vcf -f "QUAL > 20" > filteredVar.vcf

# Step 2.3 Decompose complex haplotypes
vcfallelicprimitives -k -g filteredVar.vcf > decomVar.vcf

# Step 2.4 Annotate variants
#snpEff download R64-1-1.105

snpEff ann R64-1-1.105 decomVar.vcf > annotVar.vcf
head -100 annotVar.vcf > subSampleVar.vcf

