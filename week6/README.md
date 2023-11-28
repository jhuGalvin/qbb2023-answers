# Exercise 1: Performing PCA
## Step 1.1 Compute genotype PCs
/usr/local/bin/plink_mac_20231018/plink --vcf genotypes.vcf --pca --out genotypePCs

#Exercise 2: The allele frequency spectrum
## Step 2.1 Compute allele frequencies
/usr/local/bin/plink_mac_20231018/plink --freq --vcf genotypes.vcf 

# Exercise 3: GWAS
## Step 3.1 Running the GWAS
For CB1908_IC50
/usr/local/bin/plink_mac_20231018/plink --vcf genotypes.vcf --linear --pheno CB1908_IC50.txt --covar genotypePCs.eigenvec --allow-no-sex --out CB1908_IC50_gwas_results

For GS451_IC50
/usr/local/bin/plink_mac_20231018/plink --vcf genotypes.vcf --linear --pheno GS451_IC50.txt --covar genotypePCs.eigenvec --allow-no-sex --out GS451_IC50_gwas_results

## Step 3.4 What gene could it be?
I found that this site is in the Disco-interacting protein 2 homolog B (DIP2 homolog B/DIP2B) gene. This gene may be associated with DNA methylation and AMP-binding. Where the nucleotide is present in other species, it is conserved (A, same as human). However, its Phylop conservation score is near zero (-0.556354). Further, it is within a OMIM Gene region which may be disease-causing. 