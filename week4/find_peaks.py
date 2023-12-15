#!/usr/bin/env python

import sys
from model_peaks import load_bedgraph, bin_array
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

# python find_peaks.py sample1.fwd.bg sample1.rev.bg control.fwd.bg control.rev.bg 198 sample1.peaks

def main():
    # Load file names and fragment width
    forward_fname, reverse_fname, for_ctrl_fname, rev_ctrl_fname, frag_size, output = sys.argv[1:]
    
    # Define what genomic region we want to analyze
    target = "chr2R"
    chromstart = 10000000
    chromend =  12000000
    chromlen = chromend - chromstart

    # Load the sample bedgraph data, reusing the function we already wrote
    forward = load_bedgraph(forward_fname, target, chromstart, chromend)
    reverse = load_bedgraph(reverse_fname, target, chromstart, chromend)
    
    # Load the control bedgraph data, reusing the function we already wrote
    forward_ctrl = load_bedgraph(for_ctrl_fname, target, chromstart, chromend)
    reverse_ctrl = load_bedgraph(rev_ctrl_fname, target, chromstart, chromend)

    # Combine tag densities for contol
    frag_size = int(frag_size)
    combined_ctrl = np.zeros(chromlen, float)
    combined_ctrl[frag_size//2:] += forward_ctrl[:-frag_size//2]
    combined_ctrl[:-frag_size//2] += reverse_ctrl[frag_size//2:]

    # Combine tag densities
    combined = np.zeros(chromlen, float)
    combined[frag_size//2:] += forward[:-frag_size//2]
    combined[:-frag_size//2] += reverse[frag_size//2:]

    # Adjust the control to have the same coverage as our sample
    coverage = np.sum(combined)
    ctrl_coverage = np.sum(combined_ctrl)
    ctrl_adjusted = combined_ctrl * (coverage/ctrl_coverage)

    # Create a background mean using our previous binning function and a 1K window
    # Make sure to adjust to be the mean expected per base
    binsize = 1000
    global_score = bin_array(ctrl_adjusted, binsize)/binsize

    # Find the mean tags/bp and make each background position the higher of the
    # the binned score and global background score
    mean_control_score = np.mean(ctrl_adjusted)
    background_score = np.maximum(global_score, mean_control_score)

    # Score the sample using a binsize that is twice our fragment size
    # We can reuse the binning function we already wrote
    sample_score = bin_array(combined, 2*frag_size)

    # Find the p-value for each position (you can pass a whole array of values
    # and and array of means). Use scipy.stats.poisson for the distribution.
    # Remeber that we're looking for the probability of seeing a value this large
    # or larger
    # Also, don't forget that your background is per base, while your sample is
    # per 2 * width bases. You'll need to adjust your background
    pValue_matrix =  1 - (scipy.stats.poisson.cdf(sample_score, (2*frag_size*background_score)))

    # Transform the p-values into -log10
    # You will also need to set a minimum pvalue so you doen't get a divide by
    # zero error. I suggest using 1e-250
    log_pValue_matrix = -(np.log10(pValue_matrix + 1e-250)) 

    # Write p-values to a wiggle file
    # The file should start with the line
    # "fixedStep chrom=CHROM start=CHROMSTART step=1 span=1" where CHROM and
    # CHROMSTART are filled in from your target genomic region. Then you have
    # one value per line (in this case, representing a value for each basepair).
    # Note that wiggle files start coordinates at 1, not zero, so add 1 to your
    # chromstart. Also, the file should end in the suffix ".wig"
    write_wiggle(log_pValue_matrix, target, chromstart, f"{output}.wig")

    # Write bed file with non-overlapping peaks defined by high-scoring regions 
    write_bed(log_pValue_matrix, target, chromstart, chromend, frag_size, f"{output}.bed")

# bedtools intersect -a sample1.peaks.bed -b sample2.peaks.bed > combined_peaks.bed


    # Write p-values to a wiggle file
    # The file should start with the line
    # "fixedStep chrom=CHROM start=CHROMSTART step=1 span=1" where CHROM and
    # CHROMSTART are filled in from your target genomic region. Then you have
    # one value per line (in this case, representing a value for each basepair).
    # Note that wiggle files start coordinates at 1, not zero, so add 1 to your
    # chromstart. Also, the file should end in the suffix ".wig"
def write_wiggle(pvalues, chrom, chromstart, fname):
    output = open(fname, 'w')
    print(f"fixedStep chrom={chrom} start={chromstart + 1} step=1 span=1",
          file=output)
    for i in pvalues:
        print(i, file=output)
    output.close()

# Write bed file with non-overlapping peaks defined by high-scoring regions 
def write_bed(scores, chrom, chromstart, chromend, width, fname):
    chromlen = chromend - chromstart
    output = open(fname, 'w')
    while np.amax(scores) >= 10:
        pos = np.argmax(scores)
        start = pos
        while start > 0 and scores[start - 1] >= 10:
            start -= 1
        end = pos
        while end < chromlen - 1 and scores[end + 1] >= 10:
            end += 1
        end = min(chromlen, end + width - 1)
        print(f"{chrom}\t{start + chromstart}\t{end + chromstart}", file = output)
        scores[start:end] = 0
    output.close()


if __name__ == "__main__":
    main()