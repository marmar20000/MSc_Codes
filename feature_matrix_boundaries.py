import pyBigWig
import numpy as np
import pandas as pd
#from genomicranges import GenomicRanges
from pybedtools import BedTool
import os
import sys
import ast
import multiprocessing


def create_feature_matrix(chr, plus_bw, minus_bw, features_file, window_sizes_numbers, max_feature_window,  start, end, file_histone):
    bw_plus = pyBigWig.open(plus_bw)
    bw_minus = pyBigWig.open(minus_bw)
    bw_histone = pyBigWig.open(file_histone)
    with open(features_file, "w") as features:
        while start < end:
            #print(start)
            feature_counts_pos_forward = []
            feature_counts_minus_forward = []
            feature_counts_pos_backward = []
            feature_counts_minus_backward = []
            chrom = chr
            chrom_length = bw_plus.chroms()[chrom]

            if (start + max_feature_window) < chrom_length and start - (max_feature_window) > 0:
                for repetition in window_sizes_numbers:
                    start_forwards = np.arange(start,(start + (repetition[0] * repetition[1])), repetition[1])
                    end_forwards = start_forwards + repetition[1]
                    start_backwards = np.arange((start - (repetition[0] * repetition[1])), start, repetition[1])
                    end_backwards = start_backwards + repetition[1]
                    

                    for s, e in zip(start_forwards, end_forwards):
                        pos_values = np.array(bw_plus.values(chrom, int(s), int(e)))
                        minus_values = np.array(bw_minus.values(chrom, int(s), int(e)))
                        pos_value_forward = np.abs(np.nansum(pos_values))
                        minus_value_forward = np.abs(np.nansum(minus_values))
                        feature_counts_pos_forward.append(str(pos_value_forward))
                        feature_counts_minus_forward.append(str(minus_value_forward))


                    for s, e in zip(start_backwards, end_backwards):
                        pos_values = np.array(bw_plus.values(chrom, int(s), int(e)))
                        minus_values = np.array(bw_minus.values(chrom, int(s), int(e)))
                        pos_value_backward = np.abs(np.nansum(pos_values))
                        minus_value_backward = np.abs(np.nansum(minus_values))
                        feature_counts_pos_backward.append(str(pos_value_backward))
                        feature_counts_minus_backward.append(str(minus_value_backward))


                values_of_position = np.nansum(bw_histone.values(chrom,int(start),int(int(start)+1)))
                features.write((chr)+"," + (str(start)) + "," + str(start+1)+ "," + (','.join(feature_counts_pos_forward)) + "," + (','.join(feature_counts_minus_forward)) + "," + (','.join(feature_counts_pos_backward)) + "," + (','.join(feature_counts_minus_backward)) + "," +str(values_of_position)+'\n')
            start += 10

chr = sys.argv[1]
plus_bw = sys.argv[2]
minus_bw = sys.argv[3]
features_filename = sys.argv[4]
window_sizes = ast.literal_eval(sys.argv[5])
max_feature_window = int(sys.argv[6])
start = int(sys.argv[7])
end = int(sys.argv[8])
file_histone = sys.argv[9] 
    
    

print(window_sizes, sys.argv)

create_feature_matrix(chr, plus_bw, minus_bw, features_filename, window_sizes, max_feature_window, start, end, file_histone)
