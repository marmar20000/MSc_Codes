#!/bin/bash

# This is a shared resource so please be cautious and responsible!

# 230208 edialynas v0.1 sbatch skeleton script

#SBATCH --job-name=2nd_22_10           # name of script, this is just for reporting/accounting purposes
#SBATCH --output=./2nd_22_10.out          # standard output file
#SBATCH --error=./2nd_22_10.err           # standard error file
#SBATCH --nodes=1                       # number of nodes to allocate, if your application does not run in parallel (MPI) mode set this to 1
#SBATCH --ntasks=4                     # number of cores to allocate, our nodes have 48 cores
#SBATCH --time=200:00:00                 # set a lim on the total run time, hrs:min:sec
#SBATCH --mem=200G                       # memory to allocate, our nodes have 256G, so set this to no more than 250G

# Remove the following two lines if no email notification is required
#SBATCH --mail-type=end,fail            # events to send mail on
#SBATCH --mail-user=marianna_stagaki@imbb.forth.gr    # mail recipient


WORKDIR=./                              # set this to your working directory
echo $WORKDIR
cd $WORKDIR
pwd;hostname;date                       # print-and log-working directory, hostname, start time & date, useful for job tracking and debugging

# put your pipeline commands here
python3 clean_get_informative_positions.py 'GSE89230_Normalized_PRO-seq_K562_combined_replicates_NHS_plusStrand.bigWig', 'GSE89230_Normalized_PRO-seq_K562_combined_replicates_NHS_minusStrand.bigWig', 2, 0, 100, 1000, 50, 100000,  'Info_positions_dataset1.bed'
python3 clean_get_informative_positions.py 'GSE89230_Normalized_PRO-seq_K562_combined_replicates_HS30_plusStrand_new.bigWig', 'GSE89230_Normalized_PRO-seq_K562_combined_replicates_HS30_minusStrand_new.bigWig', 2, 0, 100, 1000, 50, 100000,  'Info_positions_HS30.bed'
echo "informative positions are done"

python3 get_histone_score_remastered.py  "GSM2367733_H4ac_ChIP-seq_K562_UpState06-866_R1_36nt_NHS.bigWig" 'k562_positive_no_index.bed' 'k562_negative_no_index.bed' "Info_positions_dataset1.bed" 0.05 0.93 0.02 1_000_000 50 5 'NHS' 'True'
python3 get_histone_score_remastered.py  "GSM2367734_H4ac_ChIP-seq_K562_UpState06-866_R1_36nt_HS30_new.bigWig" 'k562_positive_no_index.bed' 'k562_negative_no_index.bed' "Info_positions_HS30.bed" 0.05 0.93 0.02 1_000_000 50 5 'HS30' 'True'

echo "histone modification scores positions are done"

python3 feature_matrix_faster.py "pro_chip_1_000_000_NHS.csv" 'GSE89230_Normalized_PRO-seq_K562_combined_replicates_NHS_plusStrand.bigWig' 'GSE89230_Normalized_PRO-seq_K562_combined_replicates_NHS_minusStrand.bigWig' "feature_matrix_1_million_NHS.bed"  '[(10,10), (10,25), (30,50), (20,500), (20,5000)]' 100000
python3 feature_matrix_faster.py "pro_chip_1_000_000_HS30.csv" 'GSE89230_Normalized_PRO-seq_K562_combined_replicates_HS30_plusStrand_new.bigWig' 'GSE89230_Normalized_PRO-seq_K562_combined_replicates_HS30_minusStrand_new.bigWig' "feature_matrix_1_million_HS30.bed"  '[(10,10), (10,25), (30,50), (20,500), (20,5000)]' 100000

echo "feature matrices is done"

python3 feature_matrix_boundaries.py chr22 GSM1480327_K562_PROseq_plus.bw GSM1480327_K562_PROseq_minus.bw feature_matrix_proseq_dataset2_chr22_per_10.bed '[(10,10), (10,25), (30,50), (20,500), (20,5000)]' 100_000 3_000_000 25_000_000 GSM4971098_mergedMasked_K4me3_0h_br1.raw.bw
