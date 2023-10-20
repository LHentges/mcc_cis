# %%
import lotron2
import csv
import pandas as pd
import numpy as np
import os
import argparse


argparser = argparse.ArgumentParser()
argparser.add_argument('--bigwig_dir', help='Directory containing bigwig files')
argparser.add_argument('--oligo_file', help='File containing oligo information - should be located in bigwig_dir')
argparser.add_argument('--write_dir', help='Directory to write output files')
argparser.add_argument('--ext', help='File extension of bigwig files - e.g. _ext_de_norm_merged.bw')
argparser.add_argument('--single_read_cov_depth', help='value of coverage for a single read after normalization', default='estimate')
argparser.add_argument('--min_read_count', help='Minimum number of reads in peak', default=5)

argparser.add_argument('--threshold_list', help='LoTron2: List of thresholds for cumulative coverage', default=[2, 4, 6])
argparser.add_argument('--background_list', help='LoTron2: List of background thresholds', default=[10**4, 10**5, 10**6])
argparser.add_argument('--window_list', help='LoTron2: List of window sizes', default=[200, 400, 600])
argparser.add_argument('--threshold_cumulative_init', help='LoTron2: Initial threshold for cumulative coverage', default=9)
argparser.add_argument('--min_size', help='LoTron2: Minimum peak width', default=50)
argparser.add_argument('--max_size', help='LoTron2: Maximum size of peaks', default=500)

args = argparser.parse_args()

bigwig_dir = args.bigwig_dir
if bigwig_dir[-1] != '/':
    bigwig_dir += '/'
oligo_file = bigwig_dir + args.oligo_file
write_dir = args.write_dir
ext = args.ext
threshold_cumulative_init = args.threshold_cumulative_init
min_size = args.min_size
max_size = args.max_size
min_read_count = args.min_read_count
single_read_cov_depth = args.single_read_cov_depth
threshold_list = args.threshold_list
background_list = args.background_list
window_list = args.window_list




# %%
def calc_peak_stats(start, end, chrom_cov_array):
    # Calculate peak stats
    peak_size = end - start
    peak_max = chrom_cov_array[start:end].max()
    peak_max_pos = chrom_cov_array[start:end].argmax() + start
    peak_sum = chrom_cov_array[start:end].sum()
    peak_mean = chrom_cov_array[start:end].mean()

    return peak_size, peak_max, int(peak_max_pos), peak_sum, peak_mean    

# %%
def estimate_single_read_cov_depth(chrom_cov_array):
    chrom_cov_nonzero = np.ma.masked_where(chrom_cov_array <= 0, chrom_cov_array)
    min_reads = chrom_cov_nonzero.min()
    return min_reads

# %%
bigwig_list = []
for file in os.listdir(bigwig_dir):
    if file.endswith(ext):
        bigwig_list.append(bigwig_dir+file)

# %%
completed_list = []
for file in os.listdir(write_dir):
    if file.endswith('.csv'):
        completed_list.append(file.split('_')[0])

# %%
exp_dict = {}
with open(oligo_file) as f:
    csv_reader = csv.reader(f, delimiter='\t')
    for row in csv_reader:
        file = bigwig_dir + row[3] + ext
        if file in bigwig_list and row[3] not in completed_list:
            if row[3] not in exp_dict:
                exp_dict[row[3]] = {}
                exp_dict[row[3]]['chr'] = row[0]
                exp_dict[row[3]]['file'] = file


print('Number of experiments found in oligo file:', len(exp_dict))
print('Number of bigwig files found:', len(bigwig_list))
print('Number of experiments with completed peak calls found:', len(completed_list))

# %%
for i, exp in enumerate(exp_dict):
    bw_lotron2 = lotron2.BigwigData(exp_dict[exp]['file'])
    chrom_cov_array = bw_lotron2.get_chrom_info_make_coverage_map(exp_dict[exp]['chr'])
    if single_read_cov_depth == 'estimate':
        single_read_cov_depth_value = estimate_single_read_cov_depth(chrom_cov_array)
    else:
        single_read_cov_depth_value = single_read_cov_depth
    min_cov = single_read_cov_depth_value * min_read_count
    mcc_peaks = bw_lotron2.call_candidate_peaks_lotron_chrom(threshold_cumulative_init, exp_dict[exp]['chr'], background_list, window_list, threshold_list, min_size, max_size, background_global_min=min_cov)
    if len(mcc_peaks) == 0:
        print(exp, 'no peaks')
        continue
    mcc_peaks[['peak_size', 'peak_max', 'peak_max_pos', 'peak_sum', 'peak_mean']] = mcc_peaks.apply(lambda row: pd.Series(calc_peak_stats(row['Start'], row['End'], chrom_cov_array)), axis=1)
    mcc_peaks['peak_max_pos'] = mcc_peaks['peak_max_pos'].astype(int)
    mcc_peaks.sort_values(by=['peak_max'], inplace=True, ascending=False)
    mcc_peaks.to_csv(write_dir + exp + '_mcc_peaks.bed', index=False, header=True, sep='\t')



