import numpy as np
import pandas as pd
import sys
import pyBigWig

def get_informative_positions(bw_file_plus, bw_file_minus,  depth_OR, depth_AND,  window_OR, window_AND, step, max_feature_window, output_file):
   
    # Open the BigWig file
    bw_plus = pyBigWig.open(bw_file_plus)
    bw_minus = pyBigWig.open(bw_file_minus)
    
    print(bw_plus.chroms())

    coordinates_above_threshold = []

    for chrom in bw_plus.chroms():
            if chrom not in  ['dm3chr2L', 'dm3chr2LHet', 'dm3chr2R','chrM',  'dm3chr2RHet', 'dm3chr3L', 'dm3chr3LHet', 'dm3chr3R', 'dm3chr3RHet',  'chrY']:
                chrom_length = bw_plus.chroms()[chrom]
                print(chrom_length)

                # Generate all start positions for the OR window (100bp)
                # In this case only 0,50 are needed as start positions since starting at 100 is a window shift ahead
                x =  np.arange(0, 100, step)
                print(x)
                for step_OR in x:
                    start_positions = np.arange((max_feature_window+step_OR-50), (chrom_length-max_feature_window-50+1), window_OR)
                    print(start_positions)
                    for start in start_positions:
                        end = start + window_OR
                        pos_value = np.array(bw_plus.values(chrom, int(start), int(end)+1))      
                        minus_value = np.array(bw_minus.values(chrom, int(start), int(end)+1))
                        if np.abs(np.nansum(pos_value)) > depth_OR or np.abs(np.nansum(minus_value)) > depth_OR:
                            coordinates_above_threshold.append((chrom, int(start+window_OR/2),  int(start+(window_OR/2)+1)))

                # Generate all start positions for the AND window (1000bp)
                # maybe we should change the step of arrange for the 1000 bp window?
                y =  np.arange(0,401, step)  
                print(y)    
                for step_AND in y:
                    start_positions = np.arange((max_feature_window+step_AND-500), (chrom_length-max_feature_window-500+1), window_AND)
                    print(start_positions)
                    for start in start_positions:
                        end = start + window_AND
                        pos_value = np.array(bw_plus.values(chrom, int(start), int(end)+1))
                        minus_value = np.array(bw_minus.values(chrom, int(start), int(end)+1))    
                        if np.abs(np.nansum(pos_value)) > depth_AND and np.abs(np.nansum(minus_value)) > depth_AND:
                            coordinates_above_threshold.append((chrom, int(start+window_AND/2),  int(start+(window_AND/2)+1)))     
                       
    bw_plus.close()
    bw_minus.close()
    df = pd.DataFrame(list(set(coordinates_above_threshold)))
    df.to_csv(output_file,sep = '\t',  index=False, header=False)

    return 'Identification of informative positions complete.'


bigwig_file_plus = sys.argv[1].split(',')[0]
bigwig_file_minus = sys.argv[2].split(',')[0]
depth_OR  =  int(sys.argv[3].split(',')[0])
depth_AND = int(sys.argv[4].split(',')[0])
window_OR = int(sys.argv[5].split(',')[0])
window_AND = int(sys.argv[6].split(',')[0])
step = int(sys.argv[7].split(',')[0])
max_feature_window = int(sys.argv[8].split(',')[0])
output_file = sys.argv[9].split(',')[0]

get_informative_positions(bigwig_file_plus, bigwig_file_minus,  depth_OR, depth_AND,  window_OR, window_AND, step, max_feature_window, output_file)