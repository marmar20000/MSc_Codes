import pyBigWig
import numpy as np
import pandas as pd
#from genomicranges import GenomicRanges
from pybedtools import BedTool
import os
import sys
import ast







def create_feature_matrix(bed_file_informative_positions, plus_bw, minus_bw, features_file, window_sizes_numbers, max_feature_window):
    bw_plus = pyBigWig.open(plus_bw)
    bw_minus = pyBigWig.open(minus_bw)
    with open(bed_file_informative_positions, "r") as bedA:
        with open(features_file, "w") as features:
            c = 0
            for chr_coordinate in bedA.readlines()[1:]:
                c += 1
                feature_counts_pos_forward = []
                feature_counts_minus_forward = []
                feature_counts_pos_backward = []
                feature_counts_minus_backward = []
                chrom = chr_coordinate.split(',')[0]
                chrom_length = bw_plus.chroms()[chrom]
                start = int(chr_coordinate.split(',')[1])
                end = int(chr_coordinate.split(',')[2])
                histone_score = (chr_coordinate.split(',')[3])

                if (start + max_feature_window) < chrom_length and start - (max_feature_window) > 0:
                    for repetition in window_sizes_numbers:
                        start_forwards = np.arange(start,(start + (repetition[0] * repetition[1])), repetition[1])
                        end_forwards = start_forwards + repetition[1]
                        start_backwards = np.arange((start - (repetition[0] * repetition[1])), start, repetition[1])
                        end_backwards = start_backwards + repetition[1]
                        
                        print("Start forwards:", start_forwards)
                        print("End forwards:", end_forwards)
                        print("Start backwards:", start_backwards)
                        print("End backwards:", end_backwards)

                        for s, e in zip(start_forwards, end_forwards):
                            pos_values = np.array(bw_plus.values(chrom, int(s), int(e)))
                            minus_values = np.array(bw_minus.values(chrom, int(s), int(e)))
                            print("Forward pos values:", pos_values)
                            print("Forward minus values:", minus_values)
                            pos_value_forward = np.abs(np.nansum(pos_values))
                            minus_value_forward = np.abs(np.nansum(minus_values))
                            feature_counts_pos_forward.append(str(pos_value_forward))
                            feature_counts_minus_forward.append(str(minus_value_forward))
                            print("Length feature_counts_pos_forward:", len(feature_counts_pos_forward))
                            print("Length feature_counts_minus_forward:", len(feature_counts_minus_forward))

                        for s, e in zip(start_backwards, end_backwards):
                            pos_values = np.array(bw_plus.values(chrom, int(s), int(e)))
                            minus_values = np.array(bw_minus.values(chrom, int(s), int(e)))
                            print("Backward pos values:", pos_values)
                            print("Backward minus values:", minus_values)
                            pos_value_backward = np.abs(np.nansum(pos_values))
                            minus_value_backward = np.abs(np.nansum(minus_values))
                            feature_counts_pos_backward.append(str(pos_value_backward))
                            feature_counts_minus_backward.append(str(minus_value_backward))
                            print("Length feature_counts_pos_backward:", len(feature_counts_pos_backward))
                            print("Length feature_counts_minus_backward:", len(feature_counts_minus_backward))

                if c % 100 == 0:
                    print(c)
                features.write((','.join(chr_coordinate.split(',')[:3])).replace(',', '_') + "," + (','.join(feature_counts_pos_forward)) + "," + (','.join(feature_counts_minus_forward)) + "," + (','.join(feature_counts_pos_backward)) + "," + (','.join(feature_counts_minus_backward)) + "," + histone_score)


pro_chip_file = sys.argv[1]
plus_bw = sys.argv[2]
minus_bw = sys.argv[3]
features_filename =  sys.argv[4]
window_sizes = ast.literal_eval(sys.argv[5])
max_feature_window =  int(sys.argv[6])
print(window_sizes, sys.argv)

create_feature_matrix(pro_chip_file, plus_bw, minus_bw, features_filename, window_sizes, max_feature_window)
