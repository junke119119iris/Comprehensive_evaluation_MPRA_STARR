# -------------------- Standard Library Imports --------------------
import argparse
import glob
import os
import sys
from multiprocessing import Pool, cpu_count
from subprocess import PIPE, Popen, STDOUT, call, run

# -------------------- Scientific Computing ------------------------
import numpy as np
import pandas as pd
import scipy

# -------------------- Bioinformatics Libraries --------------------
import pybedtools
import pyBigWig
import pysam

# -------------------- PyBedTools Temp Dir Setup -------------------
# Set temporary directory for PyBedTools operations
pybedtools.helpers.set_tempdir('/fs/cbsuhy02/storage/jz855/tmp/')

# -------------------- Project Root -------------------
# Get the current working directory as the root for relative paths
project_root = os.getcwd()

def compute_normalized_GROcap_signal(region_info):
    """
    Computes the GRO-cap signal (forward, reverse, and total) for a given genomic region.
    Signals are normalized by the region size (length in bp).

    Parameters:
    -----------
    region_info : tuple
        A tuple containing (chrom, start, end, name) representing a genomic region.

    Returns:
    --------
    list
        A list containing:
        [name,
         raw_forward_count, raw_reverse_count, raw_total_count,
         normalized_forward_signal, normalized_reverse_signal, normalized_total_signal]
    """

    chrom, start, end, name = region_info

    # Paths to forward and reverse GRO-cap bigWig signal files
    gro_cap_forward_path = os.path.join(project_root, 'data', 'Reference', 'K562_GRO_cap', 'K562_GROcap_hg38_aligned_pl.bw')
    gro_cap_forward_path = os.path.join(project_root, 'data', 'Reference', 'K562_GRO_cap', 'K562_GROcap_hg38_aligned_mn.bw')

    # Load bigWig files
    gro_cap_bw_forward = pyBigWig.open(gro_cap_forward_path)
    gro_cap_bw_reverse = pyBigWig.open(gro_cap_reverse_path)

    # Sum of GRO-cap signals across the region
    forward_count = np.nansum(gro_cap_bw_forward.values(chrom, start, end))
    reverse_count = -np.nansum(gro_cap_bw_reverse.values(chrom, start, end))  # reverse is negative strand
    total_count = forward_count + reverse_count

    # Compute region length
    size = end - start

    # Normalize signals by region size
    forward_signal = forward_count / size
    reverse_signal = reverse_count / size
    total_signal = total_count / size

    return [
        name,
        forward_count, reverse_count, total_count,
        forward_signal, reverse_signal, total_signal
    ]

def classify_region_by_grocap_signal(args):
    """
    Assigns a transcriptional activity class to a genomic region based on GRO-cap signal.

    Parameters:
    -----------
    args : tuple
        A tuple containing:
        - signal (float): normalized GRO-cap signal
        - existing_class (str): existing annotation (e.g., '.', 'bidirectional')
        - signal_only (bool): whether to assign class solely based on signal
        - low (float): lower threshold for classification
        - high (float): upper threshold for classification

    Returns:
    --------
    str
        The assigned class: 'low_transcription', 'medium_transcription', or 'high_transcription'
    """
    
    signal, existing_class, signal_only, low, high = args

    # If using signal only, ignore any existing class label
    if signal_only:
        if signal <= low:
            return 'low_transcription'
        elif low < signal <= high:
            return 'medium_transcription'
        else:
            return 'high_transcription'

    # If not using signal-only, update only when class is unassigned ('.')
    if existing_class == '.':
        if signal <= low:
            return 'low_transcription'
        elif low < signal <= high:
            return 'medium_transcription'
        else:
            return 'high_transcription'

    # Otherwise, keep the existing class
    return existing_class


