# === Standard Library ===
import os
from subprocess import call, PIPE, run
from multiprocessing import Pool, cpu_count

# === Third-Party Libraries ===
import pandas as pd
import numpy as np
import pybedtools
import pysam

# === Project Utilities ===
from utils import *

            
class BamOutToCountMatrix():
    
    """
    Combines fragment-level count data across DNA and RNA libraries
    into a unified count matrix per chromosome.
    """
    
    def __init__(self, dna_file_list, rna_file_list):
        self.dna_file_list = dna_file_list
        self.rna_file_list = rna_file_list
        self.num_rep_dna = len(self.dna_file_list)
        self.num_rep_rna = len(self.rna_file_list)
        self.sample_list = ['DNA{0}'.format(i) for i in range(1, self.num_rep_dna+1)] + ['RNA{0}'.format(i) for i in range(1, self.num_rep_rna+1)]
        
    def get_chr_list(self, chrom_size_dir):
        
        """
        Parse chromosome names from chrom.sizes file, group unplaced contigs at the end.

        Args:
            chrom_size_dir (str): Path to chrom.sizes file

        Returns:
            List[str]: Ordered chromosome names
        """
        
        chr_ref = [x.split('\t')[0] for x in open(chrom_size_dir, 'r').readlines()]
        all_other_chr = [x for x in chr_ref if '_' in x]
        chr_list = [[x] for x in chr_ref if x not in all_other_chr]
        chr_list = chr_list + [all_other_chr]
        chr_list = [y for x in chr_list for y in x]
        
        del all_other_chr
        return(chr_list)

    
    def get_count_per_sample(self, file_path, chrom, sample):
        """
        Load count entries for a specific chromosome and sample.

        Args:
            file_path (str): Path to input count BED file
            chrom (str): Chromosome name
            sample (str): Sample name to assign as column header

        Returns:
            pd.DataFrame: Indexed count table (one column)
        """
        num_row = get_row_number(file_path)
        input_file = open(file_path, 'r')

        content = []
        for r in range(0, num_row):
            line = input_file.readline().strip().split('\t')
            if(line[0] == chrom):
                content.append(line)
        input_file.close()
        
        if(len(content) > 0):
            count_mat = pd.DataFrame(content)
        else:
            count_mat = pd.DataFrame(columns=[0,1,2,3,4])

        count_mat = count_mat.astype({1:float, 2:float})
        count_mat = count_mat.astype({1:int, 2:int})

        count_mat = count_mat.set_index([0,1,2,3])
        count_mat.columns = [sample]
        
        del content
        return(count_mat)
    
    def get_count_mat_for_chr(self, args):
        """
        Aggregate counts across all samples for a given chromosome.

        Args:
            args (List[str]): [out_path, chrom]
        """
        out_path, chrom = args
        
        file_list = self.dna_file_list + self.rna_file_list
        
        for idx in range(0, len(file_list)):
            
            sample = self.sample_list[idx]
            if(idx == 0):
                count_mat = self.get_count_per_sample(file_list[idx], chrom, sample)
            else:
                curr = self.get_count_per_sample(file_list[idx], chrom, sample)
                count_mat = count_mat.join(curr, how='outer')
        
        count_mat = count_mat[self.sample_list]
        count_mat = count_mat.fillna(0)
        count_mat = count_mat.reset_index()
        
        count_mat = count_mat.astype({1:float, 2:float})
        count_mat = count_mat.astype({1:int, 2:int})
        
        if(len(count_mat) > 0):
            count_mat.to_csv(out_path+'count_mat.bed', sep='\t', mode='a', index=False, header=False)
          
        del curr, count_mat
    
    def get_count_mat(self, chrom_size_dir, out_path):
        
        """
        Process all chromosomes in parallel and build the full count matrix.

        Args:
            chrom_size_dir (str): Path to chrom.sizes file
            out_path (str): Directory to write output count matrix
        """
        # Create out_path if doesn't exist
        set_dir(out_path)
        
        chr_list = self.get_chr_list(chrom_size_dir)
        
        arg_list = [[out_path, chrom] for chrom in chr_list]
        
        with Pool(15) as pool:
            pool.map(self.get_count_mat_for_chr, arg_list)

    
class BamOutToCountMatrixGenomicBin():
    
    """
    Class to concatenate genomic bin count files from multiple DNA and RNA replicates
    into a single count matrix across all bins and samples.
    """
    
    def __init__(self, dna_file_list, rna_file_list):
        """
        Args:
            dna_file_list (List[str]): List of file paths to DNA replicate count BEDs
            rna_file_list (List[str]): List of file paths to RNA replicate count BEDs
        """
        self.dna_file_list = dna_file_list
        self.rna_file_list = rna_file_list
        self.num_rep_dna = len(self.dna_file_list)
        self.num_rep_rna = len(self.rna_file_list)
        self.sample_list = ['DNA{0}'.format(i) for i in range(1, self.num_rep_dna+1)] + ['RNA{0}'.format(i) for i in range(1, self.num_rep_rna+1)]
    
    def get_count_mat(self, out_path):
        
        """
        Merge genomic bin count matrices from DNA and RNA files into a single matrix.

        Assumes each file:
            - Has the same number of rows
            - Rows are in the same order (corresponding bins)

        Output:
            A single BED-like file: first few columns from the first file,
            followed by one column per sample (DNA and RNA), in order.
        """
        
        set_dir(out_path)
        file_list = self.dna_file_list + self.rna_file_list
        
        num_rows = get_row_number(file_list[0])
        
        open_files = [open(f, 'r') for f in file_list]
        out_file = open(out_path+'count_mat_all_bins.bed', 'w')
        
        for r in range(0, num_rows):
            for idx in range(0, len(open_files)):
                f = open_files[idx]
                if(idx == 0):
                    tmp = f.readline().strip().split('\t')
                    out_file.write('\t'.join(tmp))
                    out_file.write('\t')
                elif(idx < len(open_files)-1):
                    out_file.write(f.readline().strip().split('\t')[-1])
                    out_file.write('\t')
                else:
                    out_file.write(f.readline().strip().split('\t')[-1])
            out_file.write('\n')
            
        out_file.close()
        for f in open_files:
            f.close()
            
            
        
        
        
 