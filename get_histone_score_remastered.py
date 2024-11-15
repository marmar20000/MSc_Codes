import pyBigWig
import pandas as pd
import os
import numpy as np
import math
import pybedtools
from pybedtools import BedTool
#from memory_profiler import profile
import sys



def center_bed(input_bed_file):
    # Load the BED file
    bed = pybedtools.BedTool(input_bed_file)

    # Calculate the new start and stop positions for each interval
    def calculate_center(interval):
        center = (interval.start + interval.stop) // 2
        return center

    # Modify the intervals
    modified_intervals = []
    for interval in bed:
        center = calculate_center(interval)
        new_start = center
        new_stop = center + 1
        modified_intervals.append((interval.chrom, new_start, new_stop, interval.name, interval.score, interval.strand))

    # Create a BedTool object from the modified intervals
    return(pybedtools.BedTool(modified_intervals))

#@profile
def get_histone_peak(file_histone, bed_positive_gro_cap,  bed_negative_gro_cap, info_bed, rate_positive_positive, rate_positive_negative, rate_negative_negative, samples, step, n, condition, base_pair=True ):
    # Load the BED files and perfom merging
    info_bed = pybedtools.BedTool(info_bed)
    bed_positive_gro_cap = pybedtools.BedTool(bed_positive_gro_cap)
    bed_negative_gro_cap = pybedtools.BedTool(bed_negative_gro_cap)

    intersection_info_positive = info_bed.intersect(bed_positive_gro_cap, wa=True)
    intersection_info_negative = info_bed.intersect(bed_negative_gro_cap, wa=True)
    only_negative = bed_negative_gro_cap.subtract(info_bed)

    samples_pos_pos = np.array(intersection_info_positive.sample(n=(int(samples)*float(rate_positive_positive))))
    samples_pos_neg = np.array(intersection_info_negative.sample(n=(int(samples)*float(rate_positive_negative))))
    samples_neg_neg = only_negative.sample(n=(int(samples)*float(rate_negative_negative)))
    samples_neg_neg_base_res =  np.array(center_bed(samples_neg_neg))
    # Load the bigWig file containing histone data
    with open('pro_chip_'+samples+'_'+condition+'.csv', "w") as results_file:
        bw_histone = pyBigWig.open(file_histone)
        # Query the histone data for each peak region
        for row in range(len(samples_pos_pos)):
            if row%100==0:
                print(row)
            chrom = str(samples_pos_pos[row][0])
            start = int(int(samples_pos_pos[row][1])-int(int(step)/2))
            end = int(int(samples_pos_pos[row][1])+int(int(step)/2))
            values_of_all_region = bw_histone.values(chrom, start , end)
            bin_values = [np.nansum(values_of_all_region[i * int(n):(i + 1) * int(n)])/5 for i in range( int(step) // int(n) )]
            string_numbers = [str(num) for num in bin_values]
            result_string = ' '.join(string_numbers)
            
            info_pos_value  = bw_histone.values(str(samples_pos_pos[row][0]),int(samples_pos_pos[row][1]), int(samples_pos_pos[row][2]))
            

            if base_pair==False:
                results_file.write(str(samples_pos_pos[row][0])+","+str(samples_pos_pos[row][1])+","+str(samples_pos_pos[row][2])+","+result_string+'\n')
            else:
                results_file.write(str(samples_pos_pos[row][0])+","+str(samples_pos_pos[row][1])+","+str(samples_pos_pos[row][2])+","+str(info_pos_value)+'\n')

        for row in range(len(samples_pos_neg)):
            if row%100==0:
                print(row)
            values_of_all_region = bw_histone.values(str(samples_pos_neg[row][0]),int(int(samples_pos_neg[row][1])-int(int(step)/2)), int(int(samples_pos_neg[row][1])+int(int(step)/2)))
            bin_values = [np.nansum(values_of_all_region[i * int(n):(i + 1) * int(n)])/int(n) for i in range( int(step) // int(n)  )]
            info_pos_value  = bw_histone.values(str(samples_pos_neg[row][0]),int(samples_pos_neg[row][1]), int(samples_pos_neg[row][2]))
            string_numbers = [str(num) for num in bin_values]
            result_string = ' '.join(string_numbers)
            if base_pair==False:
                results_file.write(str(samples_pos_neg[row][0])+","+str(samples_pos_neg[row][1])+","+str(samples_pos_neg[row][2])+","+result_string+'\n')
            else:
                results_file.write(str(samples_pos_neg[row][0])+","+str(samples_pos_neg[row][1])+","+str(samples_pos_neg[row][2])+","+str(info_pos_value)+'\n')

        for row in range(len(samples_neg_neg_base_res)):
            if row%100==0:
                print(row)
            #if (int(int(samples_neg_neg_base_res[row][1])-int(int(step)/2))) >= 0 and (int(int(samples_neg_neg_base_res[row][1])+int(int(step)/2)) <= (bw_histone.chroms()[str(samples_neg_neg_base_res[row][0])])):
            values_of_all_region = bw_histone.values(str(samples_neg_neg_base_res[row][0]),int(int(samples_neg_neg_base_res[row][1])-int(int(step)/2)), int(int(samples_neg_neg_base_res[row][1])+int(int(step)/2)))
            info_pos_value  = bw_histone.values(str(samples_neg_neg_base_res[row][0]),int(samples_neg_neg_base_res[row][1]), int(samples_neg_neg_base_res[row][2]))
            #print(str(samples_neg_neg_base_res[row][0]),int(int(samples_neg_neg_base_res[row][1])-int(int(step)/2)), int(int(samples_neg_neg_base_res[row][1])+int(int(step)/2)))
            bin_values = [np.nansum(values_of_all_region[i * int(n):(i + 1) * int(n)])/int(n) for i in range( int(step) // int(n) )]
            string_numbers = [str(num) for num in bin_values]
            result_string = ' '.join(string_numbers)
            if base_pair==False:
                results_file.write(str(samples_neg_neg_base_res[row][0])+","+str(samples_neg_neg_base_res[row][1])+","+str(samples_neg_neg_base_res[row][2])+","+result_string+'\n')
            else:
                results_file.write(str(samples_neg_neg_base_res[row][0])+","+str(samples_neg_neg_base_res[row][1])+","+str(samples_neg_neg_base_res[row][2])+","+str(info_pos_value)+'\n')


        # Unload the bigWig file
        bw_histone.close()
    return "Done!"

chip_file = sys.argv[1].split(',')[0]
gro_cap_positive = sys.argv[2].split(',')[0]
gro_cap_negative =  sys.argv[3].split(',')[0]
informative_positions = sys.argv[4].split(',')[0]
rate_positive_positive=sys.argv[5].split(',')[0]
rate_positive_negative=sys.argv[6].split(',')[0]
rate_negative_negative=sys.argv[7].split(',')[0]
samples=sys.argv[8].split(',')[0]
step = sys.argv[9].split(',')[0]
n = sys.argv[10].split(',')[0]
condition =  sys.argv[11].split(',')[0]
base_pair = sys.argv[12].split(',')[0]
print(get_histone_peak(chip_file,gro_cap_positive,gro_cap_negative, informative_positions, rate_positive_positive, rate_positive_negative ,  rate_negative_negative, samples, step, n, condition, base_pair))
