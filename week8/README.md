# Exercise 1: Running Chicago

## Step 1.1 Running the runChicago.R script
Rscript runChicago.R raw/PCHIC_Data/GM_rep1.chinput,raw/PCHIC_Data/GM_rep2.chinput,raw/PCHIC_Data/GM_rep3.chinput ./chicago_Results -d ./raw/Design/ --en-feat-list ./raw/Features/featuresGM.txt -e washU_text

## Step 1.2 Feature Enrichment
CTCF - The CTCF value is about three times higher in our sample relative to the expected amount (SDOLwithSample). I am not sure if having more CTCF values means the selected region contains more TADs and is more heavily stratified/spatially controlled, or if this is indicative of some other process.

H3K4me1 - This feature is associated with gene enhancers. In our sample, this value is significantly higher than the expected amount due to random chance, indicating there is likely high enhancer activity driving interactions in the sample.

H3K4me3 & H3K27ac - These features are associated with active gene transcription and in the case of the latter, gene enhancers. These values are about 3 times higher than expected due to random chance, suggesting there is a lot of transcription happening in the sample as well, promoted by enhancers.

H3K27me3 & H3k9me3 - These features are associated with heterochromatin (gene repression). H3K27me3 levels are bout the same as expected due to random chance, while H3k9me3 levels are twice as high. Compared to enhancer and transcritpional activity, there is less repression in the sample.

## Step 2.1 Command to convert washU to UCSC file 
./convert_washU_to_UCSC.py ./raw/Design/h19_chr20and21.baitmap ./chicago_Results/data/chicago_Results_washU_text.txt chicago_Results_UCSC.bed


# Exercise 2: Visualizing interactions on the UCSC browser

## Step 2.2 Finding the most significant interactions
Here are the six top promoter-promoter interactions:
chr20 44438565 44565593 . 1390 34.77 . 0.0 chr20 44562442 44565593 PCIF1 + chr20 44438565.0 44442365.0 UBE2C +
chr20 44438565 44607204 . 1371 34.29 . 0.0 chr20 44596299 44607204 FTLP1;ZNF335 + chr20 44438565.0 44442365.0 UBE2C +
chr21 26837918 26939577 . 1360 34.02 . 0.0 chr21 26837918 26842640 snoU13 + chr21 26926437.0 26939577.0 MIR155HG +
chr20 44452862 44565593 . 1355 33.89 . 0.0 chr20 44562442 44565593 PCIF1 + chr20 44452862.0 44471524.0 SNX21;TNNC2 +
chr20 17660712 17951709 . 1354 33.85 . 0.0 chr20 17946510 17951709 MGME1;SNX5 + chr20 17660712 17672229 RRBP1 +
chr20 24972345 25043735 . 1353 33.84 . 0.0 chr20 24972345 24985047 APMAP + chr20 25036380 25043735 ACSS1 +

Here are the six top promoter-enhancer interactions:
chr21 26797667 26939577 . 1325 33.13 . 0.0 chr21 26926437 26939577 MIR155HG + chr21 26797667 26799364 . -
chr20 55957140 56074932 . 1291 32.29 . 0.0 chr20 55957140 55973022 RBM38;RP4-800J21.3 + chr20 56067414 56074932 . -
chr21 26790966 26939577 . 1166 29.17 . 0.0 chr21 26926437 26939577 MIR155HG + chr21 26790966 26793953 . -
chr20 5585992 5628028 . 1155 28.88 . 0.0 chr20 5585992 5601172 GPCPD1 + chr20 5625693 5628028 . -
chr21 26793954 26939577 . 1049 26.23 . 0.0 chr21 26926437 26939577 MIR155HG + chr21 26793954 26795680 . -
chr20 5515866 5933156 . 1043 26.08 . 0.0 chr20 5929472 5933156 MCM8;TRMT6 + chr20 5515866 5523933 . -

## Step 2.3 Visualizing top interactions
chr21:26797667-26939577; MIR155HG (MIR155HG_chr21_26797667_26939577.png) - This interaction makes sense because there is H3K4me1 at the target, meaning there may be an enhancer here. MIR155HG is responsible for regulating antiviral response and is expressed in immune cells. Its high incidence makes sense here as well, because this sample is taken from B cells.

chr20:55957140-56074932; RBM38;RP4-800J21.3 (RBM38_chr20_55957140_56074932) - This interaction is more confusing and may not represent a true promoter-enhancer pair. At the RBM30 locus, there is high transcription and all visible features, including promoter and enhancer marks in roughly the same location, which is odd. The region it is interacting with has negligible transcription and no feature annotations. There is however a nearby region of H3K4me1, meaning an enhancer could be active in the pair but that it was not correctly localized. The target, CTCFL, encodes a CTCF primarily associated with germline cells, so its pairing here is likely erroneous.

