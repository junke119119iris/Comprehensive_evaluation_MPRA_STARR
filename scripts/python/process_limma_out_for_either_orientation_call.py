# === Standard Library Imports ===
import os
import glob
import argparse
from subprocess import call, PIPE, run, Popen, STDOUT
from multiprocessing import Pool, cpu_count

# === Third-Party Library Imports ===
import pandas as pd
import numpy as np
import scipy
import pysam
import pybedtools
import statsmodels.stats.multitest as smm

# === Local Module Imports ===
from utils import *

# === Functions ===
def merge_bins(data):
    
    """
    Merge overlapping or adjacent bins from BED-like input.

    Args:
        data (pd.DataFrame): Input BED dataframe with at least ['seqnames', 'start', 'end']

    Returns:
        pd.DataFrame: Merged, sorted BED intervals
    """
    
    data = data[['seqnames', 'start', 'end']]
    data = data.sort_values(['seqnames', 'start'])
    
    bed = pybedtools.BedTool.from_dataframe(data)
    merged = bed.merge() # merge overlapping genomic regions
    merged = merged.to_dataframe(disable_auto_names=True, header=None)
    
    merged = merged.sort_values([0,1])
    
    return(merged)

def calc_activity_for_merged_peak_in_either_orientation(args):
    
    """
    Compute average activity for merged bins (regardless of strand), and extract summit info.

    Args:
        args (tuple): (idx, merged_df, raw_df)
            - idx (int): Row index of merged peak
            - merged_df (pd.DataFrame): Merged BED intervals
            - raw_df (pd.DataFrame): Raw input bin-level data

    Returns:
        list: [idx, avg_logFC, z_score_of_avg_logFC, list_of_logFCs, summit_seqname, summit_start,
               summit_end, summit_logFC, summit_z_score]
    """
    
    idx, merged, raw = args
    
    # Select the peak region
    data = merged.loc[idx, :].to_frame().transpose()
    c = data[0].tolist()[0]
    data = pybedtools.BedTool.from_dataframe(data)

    # Filter raw data for same chromosome
    raw = raw[raw['seqnames'] == c]
    raw = raw.sort_values(['seqnames', 'start'])
    raw = pybedtools.BedTool.from_dataframe(raw)

    # Intersect peak with raw bins (full overlap only)
    intersect = data.intersect(raw, wao=True, F=1)
    intersect = intersect.to_dataframe(disable_auto_names=True, header=None)
    intersect = intersect[intersect[intersect.columns.tolist()[-1]] > 0] # keep overlaps only
    intersect = intersect.drop_duplicates([6]) # remove duplicate hits
    
    # Calculate activity
    avg_activity = np.mean(intersect[7].tolist())
    avg_activity_z_score = (avg_activity-mean_neg_ctrl)/std_neg_ctrl
    activity_list = intersect[7].tolist()

    # Get the summit (bin with highest logFC)
    intersect = intersect.sort_values(7, ascending=False)
    highest_summit_idx = intersect.index.tolist()[0]

    return([idx, avg_activity, avg_activity_z_score, activity_list, 
            intersect.loc[highest_summit_idx, 3], # summit seqname
            intersect.loc[highest_summit_idx, 4], # summit start
            intersect.loc[highest_summit_idx, 5], # summit end
            intersect.loc[highest_summit_idx, 7], # summit logFC
            intersect.loc[highest_summit_idx, 12] # summit z-score
           ])


def get_activity_after_merging_in_either_orientation(merged, raw):
    
    """
    Compute merged activity and summit info for all merged peaks (either orientation).

    Args:
        merged (pd.DataFrame): Merged BED intervals
        raw (pd.DataFrame): Raw bin-level data with activity columns

    Returns:
        pd.DataFrame: Merged table with added columns:
                      mean_logFC, z_score_for_mean_logFC, logFC_list, summit info
    """
    
    raw = raw[['seqnames', 'start', 'end', 'bin_id', 'logFC', 'strand', 'P.Value', 'adj.P.Val', 'z_score', 'neg_ctrl', 'call']]
    
    raw = raw.sort_values(['seqnames', 'start'])
    
    arg_list = [[idx, merged, raw] for idx in merged.index.tolist()]
    
    with Pool(10) as pool:
        results = pool.map(calc_activity_for_merged_peak_in_either_orientation, arg_list)
        
    df = pd.DataFrame(results, columns=['idx', 'mean_logFC', 'z_score_for_mean_logFC', 'logFC_list', 
                                        'summit_seqnames', 'summit_start', 'summit_end', 'summit_logFC', 'summit_z_score'])
    df = df.set_index('idx')
    
    merged = merged.join(df, how='inner')
    
    return(merged)


def intersect_neg(query, subject, overlap_size=100):
    
    """
    Intersect two BED-like DataFrames and retain entries from `query` that 
    overlap `subject` by exactly `overlap_size` base pairs.

    Args:
        query (pd.DataFrame): BED-like DataFrame (with at least 3 columns: chrom, start, end)
        subject (pd.DataFrame): BED-like DataFrame to intersect with
        overlap_size (int): Required overlap size to retain matches (default: 100 bp)

    Returns:
        pd.DataFrame: Subset of intersected entries with exact overlap of `overlap_size`
    """
    
    # Convert both DataFrames to BedTool objects
    query = pybedtools.BedTool.from_dataframe(query)
    subject = pybedtools.BedTool.from_dataframe(subject)
    
    # Perform intersection with "write all overlaps" (wao)
    vs = query.intersect(subject, wao=True)
    # Convert back to DataFrame
    vs = vs.to_dataframe(disable_auto_names=True, header=None)
    # Filter rows with exact overlap of specified size
    vs = vs[vs[vs.columns.tolist()[-1]] == overlap_size]
    # Drop duplicate entries based on the query fragment (assumes first 5 columns describe query)
    vs = vs.drop_duplicates([0,1,2,3,4])
    
    return(vs)

# === Argument Parser ===
parser = argparse.ArgumentParser()
parser.add_argument('-p', "--working_path", help="working directory")
parser.add_argument('-o', "--out_path", help="out directory for plot")
parser.add_argument('-d', "--num_rep_DNA", type=int)
parser.add_argument('-r', "--num_rep_RNA", type=int)
parser.add_argument('-c', '--raw_count')
parser.add_argument('-t', '--tmp_dir')
parser.add_argument('--neg_ctrl_ref')

# Parse arguments
args = parser.parse_args()

# Set temporary directory for pybedtools
pybedtools.helpers.set_tempdir(args.tmp_dir)

# Set current working directory
os.chdir(args.working_path)

# Create output directory if it doesn't exist
set_dir(args.out_path)

# =============================================
# Load limma output table
limma_out = pd.read_csv('./limma_out.txt', sep='\t')
# limma_out = limma_out.set_index(['seqnames', 'start', 'end', 'strand'])
limma_out['name'] = ['bin_'+str(i) for i in range(1, len(limma_out)+1)]
limma_out['bin_id'] = ['bin_'+str(i) for i in range(1, len(limma_out)+1)]
limma_out = limma_out.set_index(['name'])

print('limma_out '+str(len(limma_out)))

# =============================================
# Annotate bins that fall into negative control reference regions

# Load negative control regions
neg_ctrl_ref = pd.read_csv(args.neg_ctrl_ref, sep='\t', header=None)
print(len(neg_ctrl_ref))

# Process forward strand bins
forward = limma_out[limma_out['strand'] == '+'][['seqnames', 'start', 'end', 'bin_id', 'strand']]
neg = neg_ctrl_ref[neg_ctrl_ref[4] == '+']
forward = intersect_neg(forward, neg)

# Process reverse strand bins
reverse = limma_out[limma_out['strand'] == '-'][['seqnames', 'start', 'end', 'bin_id', 'strand']]
neg = neg_ctrl_ref[neg_ctrl_ref[4] == '-']
reverse = intersect_neg(reverse, neg)

# Combine negative control overlaps
neg_ctrl_region = pd.concat([forward, reverse], ignore_index=True, axis=0)
print(len(neg_ctrl_region))

# Annotate limma output with negative control flag
limma_out['neg_ctrl'] = 'N'
limma_out.loc[neg_ctrl_region[3].tolist(), 'neg_ctrl'] = 'Y'

print(limma_out.groupby(['neg_ctrl']).size())

# Save negative control region activity
neg_ctrl_region = limma_out[limma_out['neg_ctrl'] == 'Y']
neg_ctrl_region.to_csv(args.out_path+'neg_ctrl_region.txt', sep='\t', index=False, header=True)
print("finished saving negative control regions")

# =============================================
# Compute Z-scores using negative control bins

# Compute mean and std from neg control logFCs
mean_neg_ctrl = np.mean(neg_ctrl_region['logFC'].values)
print(mean_neg_ctrl)

std_neg_ctrl = np.std(neg_ctrl_region['logFC'].values)
print(std_neg_ctrl)

# Create Z-score column
zscore_df = limma_out.copy()
zscore_df['z_score'] = (zscore_df['logFC']-mean_neg_ctrl)/std_neg_ctrl

# Set threshold for logFC based on Z=1.96
logFC_cutoff = 1.96*std_neg_ctrl+mean_neg_ctrl
print(logFC_cutoff)

# =============================================
# Call enhancer or repressor based on logFC and adjusted P value
zscore_df['call'] = 'inactive'
zscore_df.loc[zscore_df[(zscore_df['logFC'] >= logFC_cutoff) & (zscore_df['adj.P.Val'] < 0.05)].index.tolist(), 'call'] = 'enhancer'
zscore_df.loc[zscore_df[(zscore_df['logFC'] <= -logFC_cutoff) & (zscore_df['adj.P.Val'] < 0.05)].index.tolist(), 'call'] = 'repressor'

print(zscore_df.groupby(['call']).size())

# Keep only selected columns and reformat index

zscore_df = zscore_df.set_index(['seqnames', 'start', 'end', 'strand'])
zscore_df = zscore_df[['logFC', 'P.Value', 'adj.P.Val', 'bin_id', 'neg_ctrl', 'z_score', 'call']]
zscore_df.head(1)

# =============================================
# If enhancer bins exist, join with raw counts and export per-bin data

if(len(zscore_df[zscore_df['call'] == 'enhancer']) > 0):
    
    # Read raw count matrix if provided
    if(args.raw_count != None):
        raw_count = pd.read_csv(args.raw_count, sep='\t', header=None)
        raw_count.columns = ['seqnames', 'start', 'end', 'strand'] + ['DNA{0}'.format(i) for i in range(1, args.num_rep_DNA+1)] + ['RNA{0}'.format(i) for i in range(1, args.num_rep_RNA+1)]
    else:
        raw_count = './filtered_raw_count.txt'
        raw_count = pd.read_csv(raw_count, sep='\t')
        raw_count.columns = ['seqnames', 'start', 'end', 'strand'] + ['DNA{0}'.format(i) for i in range(1, args.num_rep_DNA+1)] + ['RNA{0}'.format(i) for i in range(1, args.num_rep_RNA+1)]
    
    raw_count = raw_count.set_index(['seqnames', 'start', 'end', 'strand'])
    print('raw_count '+str(len(raw_count)))
    
    # Join raw counts to annotated bin-level results
    bin_level_result_all = zscore_df.join(raw_count, how='left')
    bin_level_result_all = bin_level_result_all.reset_index()
    bin_level_result_all = bin_level_result_all[['seqnames', 'start', 'end', 'bin_id', 'logFC', 'strand', 'P.Value', 'adj.P.Val', 'z_score', 'neg_ctrl', 'call']+['DNA{0}'.format(i) for i in range(1, args.num_rep_DNA+1)] + ['RNA{0}'.format(i) for i in range(1, args.num_rep_RNA+1)]]
    
    # Save bin-level results
    set_dir(args.out_path+'bin_level/')
    bin_level_result_all.to_csv(args.out_path+'bin_level/bin_level_all_result.bed.gz', sep='\t', index=False, header=False, compression='gzip')
    print('finished saving bin-level results in each orientation')
    
    # =============================================
    # Merge enhancer and repressor bins into peaks
    
    enhancer_bins = zscore_df[zscore_df['call'] == 'enhancer'].reset_index()
    repressor_bins = zscore_df[zscore_df['call'] == 'repressor'].reset_index()
    
    set_dir(args.out_path+'merged_peak/')
    
    for data, fname in zip([enhancer_bins, repressor_bins], 
                           ['merged_enhancer_peak_in_either_orientation.bed.gz', 
                            'merged_repressor_peak_in_either_orientation.bed.gz']):
        
        # Merge adjacent/overlapping bins into peaks
        merged = merge_bins(data)
        merged = get_activity_after_merging_in_either_orientation(merged, data)
        
        # Add required columns
        merged['strand'] = '.'
        merged['size'] = merged[2].values - merged[1].values
        
        merged.columns = ['seqnames', 'start', 'end', 'mean_logFC', 'z_score_for_mean_logFC', 'logFC_list', 
                          'summit_seqnames', 'summit_start', 'summit_end', 'summit_logFC', 'summit_z_score', 'strand', 'size']
        merged = merged.sort_values(['seqnames', 'start'])
        
        # Add peak names and reorder
        merged['name'] = ['peak{0}'.format(i) for i in range(1, len(merged)+1)]
        merged = merged[['seqnames', 'start', 'end', 'name', 'mean_logFC', 'strand', 
                         'z_score_for_mean_logFC', 'logFC_list', 
                         'summit_seqnames', 'summit_start', 'summit_end', 'summit_logFC', 'summit_z_score', 'size']]
        
        # Write merged peak file
        merged.to_csv(args.out_path+'merged_peak/'+fname, sep='\t', index=False, header=False, compression='gzip')
        print('finished saving ' + fname)