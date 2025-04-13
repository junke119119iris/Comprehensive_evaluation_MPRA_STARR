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

def calc_activity_for_merged_peak_both_orientation(args):
    
    """
    Calculate activity for a merged peak (bin) by intersecting it with strand-specific bins,
    using both forward and reverse orientation activity.

    Args:
        args (tuple): (idx, merged, raw)

    Returns:
        list: [
            idx, mean activity, z-score, list of activities,
            summit_chr, summit_start, summit_end, summit_logFC, summit_z_score
        ]
    """
    
    idx, merged, raw = args
    
    # Extract one merged peak
    data = merged.loc[idx, :].to_frame().transpose()
    c = data[0].tolist()[0]
    data = pybedtools.BedTool.from_dataframe(data)

    # Subset raw bins for the chromosome of interest
    raw = raw[raw['seqnames'] == c]
    raw = raw.sort_values(['seqnames', 'start'])
    raw = pybedtools.BedTool.from_dataframe(raw)

    # Intersect the merged peak with bins
    intersect = data.intersect(raw, wao=True, F=1)
    intersect = intersect.to_dataframe(disable_auto_names=True, header=None)
    
    # Filter out bins with no overlap
    intersect = intersect[intersect[intersect.columns.tolist()[-1]] > 0]
    intersect = intersect.drop_duplicates([6]) # remove duplicates
    
    # Combine forward and reverse logFC columns
    activity_list = intersect[11].tolist() + intersect[17].tolist()
    avg_activity = np.mean(activity_list)

    # Z-score using negative control distribution
    avg_activity_z_score = (avg_activity-mean_neg_ctrl)/std_neg_ctrl

    # Determine summit (bin with max forward logFC)
    intersect = intersect.sort_values(7, ascending=False)
    highest_summit_idx = intersect.index.tolist()[0]

    return([idx, avg_activity, avg_activity_z_score, activity_list, 
            intersect.loc[highest_summit_idx, 3], # summit seqname
            intersect.loc[highest_summit_idx, 4], # summit start
            intersect.loc[highest_summit_idx, 5], # summit end
            intersect.loc[highest_summit_idx, 7], # summit logFC
            intersect.loc[highest_summit_idx, 9]  # summit z-score
           ])


def get_activity_after_merging_both_orientation(merged, raw):
    
    """
    Compute activity across merged bins using both forward and reverse strand values.

    Args:
        merged (pd.DataFrame): Merged BED intervals
        raw (pd.DataFrame): BED bins with logFC/z-score per strand

    Returns:
        pd.DataFrame: Merged table with activity summary columns appended
    """
    
    # Select required columns only
    raw = raw[['seqnames', 'start', 'end', 'name', 'activity', 'strand', 
             'activity_z_score', 'call', 
             'logFC_for', 'P.Value_for', 'q_for', 'z_score_for', 'neg_ctrl_for', 'call_for',  
             'logFC_rev', 'P.Value_rev', 'q_rev', 'z_score_rev', 'neg_ctrl_rev', 'call_rev']]
    raw = raw.sort_values(['seqnames', 'start'])
    arg_list = [[idx, merged, raw] for idx in merged.index.tolist()]
    
    # Parallel activity calculation
    with Pool(10) as pool:
        results = pool.map(calc_activity_for_merged_peak_both_orientation, arg_list)
        
    # Merge results back to merged BED
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

# # Keep only selected columns and reformat index

zscore_df = zscore_df.set_index(['seqnames', 'start', 'end', 'strand'])
zscore_df = zscore_df[['logFC', 'P.Value', 'adj.P.Val', 'bin_id', 'neg_ctrl', 'z_score', 'call']]
zscore_df.head(1)

# =============================================
# If enhancer bins exist, proceed to orientation-independent analysis

if(len(zscore_df[zscore_df['call'] == 'enhancer']) > 0):

    zscore_df = zscore_df.reset_index()
    
    # Separate forward and reverse calls
    forward = zscore_df[zscore_df['strand'] == '+']
    forward = forward.set_index(['seqnames', 'start', 'end'])
    forward = forward[['logFC', 'P.Value', 'z_score', 'neg_ctrl', 'call']]
    forward.columns = [x+'_for' for x in forward.columns.tolist()]
    
    reverse = zscore_df[zscore_df['strand'] == '-']
    reverse = reverse.set_index(['seqnames', 'start', 'end'])
    reverse = reverse[['logFC', 'P.Value', 'z_score', 'neg_ctrl', 'call']]
    reverse.columns = [x+'_rev' for x in reverse.columns.tolist()]
    
    # Join both orientations and keep only complete pairs
    both = forward.join(reverse, how='outer')
    print(len(both))
    both = both.dropna()
    print(len(both))
    
    # Average activity and compute z-score
    both['activity'] = (both['logFC_for'].values + both['logFC_rev'].values)/2
    both['activity_z_score'] = (both['activity']-mean_neg_ctrl)/std_neg_ctrl
    
    # Set up bin name and strand
    both = both.reset_index()
    both = both.sort_values(['seqnames', 'start'])
    both['name'] = ['bin_tested_both_{0}'.format(i) for i in range(1, len(both)+1)]
    both['strand'] = '.'
    both = both.set_index('name')
    
    # Recalculate FDRs for both strands
    tmp1 = both[['P.Value_for']]
    tmp1.columns = ['p']
    tmp1['type'] = 'for'
    tmp1 = tmp1.reset_index()
    
    tmp2 = both[['P.Value_rev']]
    tmp2.columns = ['p']
    tmp2['type'] = 'rev'
    tmp2 = tmp2.reset_index()
    
    tmp = pd.concat([tmp1, tmp2], axis=0)
    tmp['adj.p'] = smm.fdrcorrection(tmp['p'].values, alpha=0.05, method='indep', is_sorted=False)[1]

    # Split corrected p-values
    tmp1 = tmp[tmp['type'] == 'for']
    tmp1 = tmp1.set_index('name')
    tmp1 = tmp1[['adj.p']]
    tmp1.columns = ['q_for']
    
    tmp2 = tmp[tmp['type'] == 'rev']
    tmp2 = tmp2.set_index('name')
    tmp2 = tmp2[['adj.p']]
    tmp2.columns = ['q_rev']

    # Join back to main
    both = both.join(tmp1, how='inner')
    print(len(both))
    both = both.join(tmp2, how='inner')
    print(len(both))
    
    both = both.reset_index()
    
    # Final enhancer/repressor calls (requires both strands to be significant)
    both['call'] = 'inactive'
    both.loc[both[(both['logFC_for'] >= logFC_cutoff) & (both['q_for'] < 0.05) & 
                  (both['logFC_rev'] >= logFC_cutoff) & (both['q_rev'] < 0.05)].index.tolist(), 'call'] = 'enhancer'
    both.loc[both[(both['logFC_for'] <= -logFC_cutoff) & (both['q_for'] < 0.05) & 
                  (both['logFC_rev'] <= -logFC_cutoff) & (both['q_rev'] < 0.05)].index.tolist(), 'call'] = 'repressor'
    print(both.groupby(['call']).size())
    
    # Final columns
    both = both[['seqnames', 'start', 'end', 'name', 'activity', 'strand', 
             'activity_z_score', 'call', 
             'logFC_for', 'P.Value_for', 'q_for', 'z_score_for', 'neg_ctrl_for', 'call_for',  
             'logFC_rev', 'P.Value_rev', 'q_rev', 'z_score_rev', 'neg_ctrl_rev', 'call_rev']]
    both.to_csv(args.out_path+'bin_level/all_bin_tested_in_both_orientations.bed.gz', sep='\t', index=False, header=False, compression='gzip')
    print('finished saving all bin-level results for bins tested in both orientations')
    
    # ===============================
    # Merge orientation-independent enhancers into peak regions
    # ===============================
    both_up = both[both['call'] == 'enhancer']
    print(len(both_up))
    both_up = both_up.sort_values(['seqnames', 'start'])
    
    merged = merge_bins(both_up)
    merged = get_activity_after_merging_both_orientation(merged, both_up)
    
    merged.columns = ['seqnames', 'start', 'end', 'mean_logFC', 'z_score_for_mean_logFC', 'logFC_list', 
                  'summit_seqnames', 'summit_start', 'summit_end', 'summit_logFC', 'summit_z_score']
    merged = merged.sort_values(['seqnames', 'start'])
    
    merged['name'] = ['peak_both_{0}'.format(i) for i in range(1, len(merged)+1)]
    merged['strand'] = '.'
    merged['size'] = merged['end'].values - merged['start'].values
    merged = merged[['seqnames', 'start', 'end', 'name', 'mean_logFC', 'strand', 'z_score_for_mean_logFC', 'logFC_list', 
                     'summit_seqnames', 'summit_start', 'summit_end', 'summit_logFC', 'summit_z_score', 'size']]
    merged.to_csv(args.out_path+'merged_peak/merged_enhancer_peak_orientation_independent.bed.gz', sep='\t', header=False, index=False, compression='gzip')
    print(len(merged))
    print('finished saving orientation-independent merged peaks')

    # ===============================
    # Repeat for repressor bins
    # ===============================
    both_down = both[both['call'] == 'repressor']
    print(len(both_down))
    both_down = both_down.sort_values(['seqnames', 'start'])
    
    # merge orientation-independent peaks
    merged = merge_bins(both_down)
    merged = get_activity_after_merging_both_orientation(merged, both_down)
    
    merged.columns = ['seqnames', 'start', 'end', 'mean_logFC', 'z_score_for_mean_logFC', 'logFC_list', 
                  'summit_seqnames', 'summit_start', 'summit_end', 'summit_logFC', 'summit_z_score']
    merged = merged.sort_values(['seqnames', 'start'])
    
    merged['name'] = ['peak_both_{0}'.format(i) for i in range(1, len(merged)+1)]
    merged['strand'] = '.'
    merged['size'] = merged['end'].values - merged['start'].values
    merged = merged[['seqnames', 'start', 'end', 'name', 'mean_logFC', 'strand', 'z_score_for_mean_logFC', 'logFC_list', 
                     'summit_seqnames', 'summit_start', 'summit_end', 'summit_logFC', 'summit_z_score', 'size']]
    merged.to_csv(args.out_path+'merged_peak/merged_repressor_peak_orientation_independent.bed.gz', sep='\t', header=False, index=False, compression='gzip')
    print(len(merged))
    print('finished saving orientation-independent merged peaks')
    
    
    
    
    



