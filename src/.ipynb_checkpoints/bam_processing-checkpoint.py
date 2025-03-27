# === Standard Library ===
import os
from subprocess import call, PIPE, run
from multiprocessing import Pool, cpu_count

# === Data Handling ===
import pandas as pd
import numpy as np

# === Bioinformatics Tools ===
import pysam
import pybedtools

# === Project Utilities ===
from utils import *

class ProcessBamSTARR:
    
    """
    Class for processing BAM files into sorted BED fragments or genomic bins suitable for downstream analysis.

    Attributes:
        bam_dir (str): Path to input BAM file.
        min_frag_size (int): minimal of fragment size to be kept for downstream analysis
        max_frag_size (int): maximum of fragment size to be kept for downstream analysis
        prefix (str): Prefix of directory to save output.
        tmp_dir (str): Path to directory to save temporary files. 
        sample (str): Sample name used for outpuot file naming. 
    """
    
    def __init__(self, bam_dir, min_frag_size, max_frag_size, prefix, tmp_dir):#, remove_duplicate, collapse):
        self.bam_dir = bam_dir
        self.min_frag_size = min_frag_size
        self.max_frag_size = max_frag_size
        self.prefix = prefix
        self.tmp_dir = tmp_dir
        self.sample = self.prefix.split('/')[-2]
        set_pybedtools_tmp_dir(self.tmp_dir)
        
    def get_chr_list(self, chrom_size_dir):
        
        """
        Get the list of chromosomes to process from a BAM file, filtered by reference.

        This method reads the list of chromosomes from a reference file (e.g., chrom.sizes),
        and intersects it with the chromosomes that actually have aligned reads in the BAM file.

        Args:
            chrom_size_dir (str): Path to a chromosome size file (e.g., UCSC chrom.sizes format)

        Returns:
            list: List of chromosomes that are present in both the reference and the BAM file,
                  and have at least one aligned read.
        """
        
        # Get reference chromosome names from chrom size file
        chr_ref = [x.split('\t')[0] for x in open(chrom_size_dir, 'r').readlines()]
        
        # Open BAM file and get statistics on each chromosome
        b = pysam.AlignmentFile(self.bam_dir, 'rb')
        chr_list = b.get_index_statistics()
        
        # Filter chromosomes with at least one mapped read
        chr_list = [chr_list[x][0] for x in range(0, len(chr_list)) if chr_list[x][1] > 0]
        b.close()
        
        # Return only chromosomes that are present in both the BAM and reference list
        return(list(set(chr_list).intersection(set(chr_ref))))
     
    def process_bam_per_chrom(self, args):
        
        """
        Process a BAM file for a given chromosome to extract valid fragment pairs,
        record their orientation, and write them to a BED file.

        This method processes paired-end reads and filters fragments based on:
        - mapping quality
        - fragment size
        - proper pairing
        - optional duplicate removal

        It separates fragments into forward and reverse orientations,
        and records statistics on the number of fragments seen, filtered, and written.

        Args:
            args (tuple): (chr, remove_duplicate)
                chr (str): Chromosome name to process
                remove_duplicate (bool): Whether to remove duplicate reads

        Returns:
            tuple:
                (frag_count_fwd, frag_count_rev,
                 frag_count_ready_for_use_fwd, frag_count_ready_for_use_rev,
                 frag_count_used_fwd, frag_count_used_rev)
        """

        chr, remove_duplicate = args
    
        # Open BAM file
        b = pysam.AlignmentFile(self.bam_dir, 'rb')

        # Fragment orientation counters 
        frag_count_fwd = 0
        frag_count_rev = 0
        frag_count_ready_for_use_fwd = 0
        frag_count_ready_for_use_rev = 0
        frag_count_used_fwd = 0
        frag_count_used_rev = 0
        
        # Output BED file for fragments
        with open(self.prefix+'alignment.bed', 'a') as s:
            r_cache = {} # Temporary storage for paired-end reads

            for read in b.fetch(reference=chr):
                rid = read.query_name
                
                # Filtering valid reads (no secondary alignments, good MAPQ, etc.)
                if(read.is_proper_pair and not read.is_secondary and read.mapping_quality >= 10 and not read.has_tag("SA")):

                    if(rid in r_cache):
                        # Get the paired reads
                        if(read.is_read2):
                            read1 = r_cache[rid]
                            read2 = read
                        else:
                            read1 = read
                            read2 = r_cache[rid]
                        del r_cache[rid]
                        
                        # Coordinates and lengths for each read
                        r1_c = b.get_reference_name(read1.reference_id)
                        r1_s = read1.reference_start
                        r1_e = read1.reference_end
                        r1_l = read1.template_length

                        r2_c = b.get_reference_name(read2.reference_id)
                        r2_s = read2.reference_start
                        r2_e = read2.reference_end
                        r2_l = read2.template_length
                        
                        # Ensure both reads map to the same chromosome
                        if(r1_c == r2_c):

                            # Forward strand: read1 is forward, read2 is reverse
                            if(not read1.is_reverse and read2.is_reverse):
                                frag_count_fwd += 1

                                # Filter based on size 
                                if(int(self.min_frag_size) <= r1_l and int(self.max_frag_size) >= r1_l):

                                    # Make sure reads are properly paired
                                    if(r2_e - r1_s == r1_l):
                                        frag_count_ready_for_use_fwd += 1
                                        
                                        # Remove duplicate if necessary
                                        if(remove_duplicate):
                                            if(read.is_duplicate):
                                                continue
                                                
                                        frag_count_used_fwd += 1
                                        # Write fragment to BED
                                        s.write("%s\t%i\t%i\t%s\t.\t%s\n" % (chr, r1_s, r2_e, rid, "+"))
                                        
                            # Reverse strand: read1 is reverse, read2 is forward              
                            elif(read1.is_reverse and not read2.is_reverse):
                                frag_count_rev += 1
                
                                # Filter based on size 
                                if(int(self.min_frag_size) <= r2_l and int(self.max_frag_size) >= r2_l):
                                    # Make sure reads are properly paired
                                    if(r1_e - r2_s == r2_l):
                                        
                                        frag_count_ready_for_use_rev += 1
                                        
                                        # Remove duplicate if necessary
                                        if(remove_duplicate):
                                            if(read.is_duplicate):
                                                continue
                                                
                                        frag_count_used_rev += 1

                                        # Write fragment to BED
                                        s.write("%s\t%i\t%i\t%s\t.\t%s\n" % (chr, r2_s, r1_e, rid, "-"))
                                        
                    else:
                        # Cache the reaad until its pair shows up
                        r_cache[rid] = read
        # Clean up
        del r_cache
        b.close()
        
        return(frag_count_fwd, frag_count_rev, 
               frag_count_ready_for_use_fwd, frag_count_ready_for_use_rev, 
               frag_count_used_fwd, frag_count_used_rev)
    
    
    def process_bam(self, chrom_size_dir, remove_duplicate):
        """
        Process BAM file to extract valid fragments from all chromosomes using parallelization.

        This function:
        - Retrieves chromosomes to process using reference + BAM content
        - Calls `process_bam_per_chrom()` in parallel for each chromosome
        - Aggregates fragment counts from both orientations
        - Sorts and compresses the final BED output

        Args:
            chrom_size_dir (str): Path to chromosome size reference file.
            remove_duplicate (bool): Whether to remove PCR duplicates.

        Returns:
            None
        """
        
        # Get list of chromosomes that exist in both reference and BAM file
        chr_list = self.get_chr_list(chrom_size_dir)
        
        # Prepare arguments for multiprocessing
        arg_list = [[chr, remove_duplicate] for chr in chr_list]
        
        # Parallel BAM processing per chromosome
        with Pool(15) as pool:
            list_counts = pool.map(self.process_bam_per_chrom, arg_list)
            
        # Aggregate fragment statistics across all chromosomes
        frag_count_fwd, frag_count_rev, frag_count_ready_for_use_fwd, frag_count_ready_for_use_rev, frag_count_used_fwd, frag_count_used_rev = np.sum(list_counts, axis=0)
        total_used_frag_count = frag_count_used_fwd + frag_count_used_rev
        
        # Print fragment statistics
        print("%s fragments are extracted from forward orientation" % (frag_count_used_fwd))
        print("%s fragments are extracted from reverse orientation" % (frag_count_used_rev))
        print("%s fragments in total are extracted" % (total_used_frag_count))
        total_ready_for_used_frag_count = frag_count_ready_for_use_fwd + frag_count_ready_for_use_rev
        print("%s fragments in total (including duplicates)" % (total_ready_for_used_frag_count))
        
        # Sort, clean and compress final alignment BED file
        safe_bedsort(self.prefix+'alignment.bed', self.prefix+'alignment.sorted.bed')
        safe_remove(self.prefix+'alignment.bed')
        gzip_file(self.prefix+'alignment.sorted.bed')

        
    def get_frag_count(self, collapse):
        
        """
        Count number of fragments from sorted alignment BED file.

        This function:
        - Reads the compressed or uncompressed alignment BED file
        - Optionally collapses all counts to 1 (for presence/absence analysis)
        - Groups by chrom/start/end/strand to count fragment support
        - Writes the resulting count matrix to `collapsed/` or `all/` directory

        Args:
            collapse (bool): If True, collapse all fragment counts to 1.

        Returns:
            None
        """

        # Determine output directory
        if(collapse):
            out_path = self.prefix+'collapsed/'
        else:
            out_path = self.prefix+'all/'
            
        # Unzip file if needed
        if(os.path.exists(self.prefix + 'alignment.sorted.bed'+'.gz') and not os.path.exists(self.prefix + 'alignment.sorted.bed')):
            gunzip_file(self.prefix + 'alignment.sorted.bed'+'.gz', self.prefix + 'alignment.sorted.bed', True)
    
        # Read large BED file in chunks to conserve memory
        chunk_list = []
        chunksize = 10**6
        for chunk in pd.read_csv(self.prefix + 'alignment.sorted.bed', sep='\t', chunksize = chunksize, header=None):
            chunk_list.append(chunk)
        data = pd.concat(chunk_list)
        del chunk_list
        
        # Group by chrom, start, end, strand and count reads
        count = data.groupby([0,1,2,5]).size().to_frame()
        count.columns = [self.sample]
        count = count.reset_index()
        
        # Collapse all counts to 1 if requested
        if(collapse):
            count[self.sample] = 1
            
        # Write count matrix to BED file
        count.to_csv(out_path+'count.bed', sep='\t', index=False, header=False)
        
        # Clean up: re-remove unzipped BED file if it was temporarily created
        if(os.path.exists(self.prefix + 'alignment.sorted.bed'+'.gz') and os.path.exists(self.prefix + 'alignment.sorted.bed')):
            safe_remove(self.prefix + 'alignment.sorted.bed')
            
        del data, count
        
    def get_binned_frag_count(self, binSize, collapse):
        
        """
        Count fragments grouped into fragment bins of a given size.

        This function:
        - Unzips sorted alignment BED if needed
        - Optionally collapses duplicate fragments
        - Assigns fragments into genomic bins based on start/end positions
        - Groups and counts fragments by binned region and strand
        - Saves result as a BED-like file for downstream usage

        Args:
            binSize (int): Size of genomic bins (e.g. 500 for 500bp bins)
            collapse (bool): If True, collapse duplicate fragments

        Returns:
            None
        """
        
        # Determine output directory
        if(collapse):
            out_path = self.prefix+'collapsed/'
        else:
            out_path = self.prefix+'all/'
        
        # Unzip alignment file if needed
        if(os.path.exists(self.prefix + 'alignment.sorted.bed'+'.gz') and not os.path.exists(self.prefix + 'alignment.sorted.bed')):
            gunzip_file(self.prefix + 'alignment.sorted.bed'+'.gz', self.prefix + 'alignment.sorted.bed', True)
        
        # Load fragments as a pandas DataFrame via pyBedTools
        align =  pybedtools.BedTool(self.prefix+'alignment.sorted.bed')
        align = align.to_dataframe(disable_auto_names=True, header=None)
        
        # Collapse to unique fragments if requested
        if(collapse):
            align = align[[0,1,2,5]].drop_duplicates()
        
        # Assign binned start/end based on bin size (rounded)
        align['binned_start'] = binSize*np.round(align[1].values/binSize)+1
        align['binned_end'] = binSize*np.round(align[2].values/binSize)+1
        
        # Groupby chrom, binned start, end, and strand; then count
        count = align.groupby([0, 'binned_start', 'binned_end', 5]).size().to_frame()
        count.columns = [self.sample]
        count = count.reset_index()
    
        # Wrtie count matrix to BED file
        count.to_csv(out_path + 'binned_frag_count_binSize_'+str(binSize)+'.bed', sep='\t', index=False, header=False)
        
        # Clean up temporary unzipped file if it was created
        if(os.path.exists(self.prefix + 'alignment.sorted.bed'+'.gz') and os.path.exists(self.prefix + 'alignment.sorted.bed')):
            safe_remove(self.prefix + 'alignment.sorted.bed')
        
        del align, count
        
        
    def make_bin(self, windowSize, stepSize, chrom_size_dir):
        
        """
        Generate sliding genomic windows using specified window and step sizes.

        This function uses pybedtools to create genome-wide bins/sliding windows
        and outputs a sorted BED file for downstream fragment binning.

        Args:
            windowSize (int): Length of each genomic window (e.g., 500 for 500bp)
            stepSize (int): Step between adjacent windows (e.g., 10 for 10bp step size)
            chrom_size_dir (str): Path to chromosome size file (UCSC format)

        Returns:
            None
        """
        
        # Determine prefix to save output one level above
        prefix = '/'.join(self.prefix.split('/')[:-3])+'/'
        
        # Create sliding windows using pybedtools
        windows = pybedtools.BedTool().window_maker(g=chrom_size_dir, w=windowSize, s=stepSize)
        windows = windows.to_dataframe(disable_auto_names=True, header=None)

        # Save raw unsorted BED file
        tmp_bed = prefix+str(windowSize)+'_'+str(stepSize)+'_tmp.windows.bed'
        windows.to_csv(tmp_bed, sep='\t', index=False, header=False)
        del windows
        
        # Sort the BED file
        sorted_bed = prefix+str(windowSize)+'_'+str(stepSize)+'_tmp.windows.sorted.bed'
        safe_bedsort(tmp_bed, sorted_bed)
        
        # Remove unsorted temp file
        safe_remove(tmp_bed)
        
    def count_for_bin(self, windows, query):
        
        """
        Count number of fragments (reads) that fully overlap each genomic window.

        This function uses pybedtools to calculate read depth by computing
        full-overlap coverage of query fragments across a set of pre-defined bins.

        Args:
            windows (pybedtools.BedTool): BED intervals representing genomic windows
            query (pybedtools.BedTool): BED intervals representing aligned fragments

        Returns:
            pandas.DataFrame: A dataframe with columns:
                - 'seqnames': Chromosome
                - 'start': Start coordinate of window
                - 'end': End coordinate of window
                - 'count': Number of reads/fragments fully overlapping the window
        """
        
        # Compute full overlap (f=1.0) coverage for each window
        depth = windows.coverage(query, sorted=True, counts=True, f=1.0)
        depth = depth.to_dataframe(disable_auto_names=True, header=None)
        
        # Rename columns for clarity
        depth.columns = ['seqnames', 'start', 'end', 'count']
        
        # Cleanup references
        del windows, query
        
        return(depth)

    def get_genomic_bin_count(self, windowSize, stepSize, collapse, chrom_size_dir):
        
        """
        Compute fragment count per genomic bin using full-overlap criteria,
        separated by strand (+/-) and optionally collapsed by fragment uniqueness.

        This function:
        - Generates sliding windows if not already present
        - Prepares strand-separated fragment BED files
        - Counts full-overlap fragments per window and strand
        - Outputs sorted genomic bin count BED file for downstream analysis

        Args:
            windowSize (int): Size of each genomic window
            stepSize (int): Step between adjacent windows
            collapse (bool): Whether to collapse duplicate fragments
            chrom_size_dir (str): Path to chromosome size file

        Returns:
            None
        """

        # Determine output directory based on collapse mode
        if(collapse):
            out_path = self.prefix+'collapsed/'
        else:
            out_path = self.prefix+'all/'
            
        # Create windows if not aready done
        if(not os.path.exists('/'.join(self.prefix.split('/')[:-3])+'/'+str(windowSize)+'_'+str(stepSize)+'_tmp.windows.sorted.bed')):
            self.make_bin(windowSize, stepSize, chrom_size_dir)
        
        # Prepare strand-separated fragment BEDs if not already present
        forward_file = out_path+'alignment.forward.sorted.bed'
        reverse_file = out_path+'alignment.reverse.sorted.bed'  
        
        if(not os.path.exists(forward_file) and not os.path.exists(reverse_file)):

            # Unzip alignment if needed
            if(os.path.exists(self.prefix + 'alignment.sorted.bed'+'.gz') and not os.path.exists(self.prefix + 'alignment.sorted.bed')):
                gunzip_file(self.prefix + 'alignment.sorted.bed'+'.gz', self.prefix + 'alignment.sorted.bed', True)

            # Read alignment and filter necessary columns
            align =  pybedtools.BedTool(self.prefix+'alignment.sorted.bed').to_dataframe(disable_auto_names=True, header=None)
            if(collapse):
                align = align[[0,1,2,5]].drop_duplicates()
            else:
                align = align[[0,1,2,5]]
            align.columns = [0,1,2,3] # Rename for pybedtools compatibility
            
            # Split into forward and reverse strand
            forward = align[align[3] == '+']
            reverse = align[align[3] == '-']
            
            # Write to temporary strand-separated files
            forward.to_csv(forward_file, sep='\t', index=False, header=False)
            reverse.to_csv(reverse_file, sep='\t', index=False, header=False)

            del forward, reverse, align
  
        # Load input files
        windows = pybedtools.BedTool('/'.join(self.prefix.split('/')[:-3])+'/'+str(windowSize)+'_'+str(stepSize)+'_tmp.windows.sorted.bed')
        forward = pybedtools.BedTool(forward_file)
        reverse = pybedtools.BedTool(reverse_file)
        
        # Count full-overlap coverage by strand
        forward_depth = self.count_for_bin(windows, forward, full=True)
        reverse_depth = self.count_for_bin(windows, reverse, full=True)
            
        forward_depth['strand'] = '+'
        reverse_depth['strand'] = '-'
        
        # Reorder and merge
        forward_depth = forward_depth[['seqnames', 'start', 'end', 'strand', 'count']]
        reverse_depth = reverse_depth[['seqnames', 'start', 'end', 'strand', 'count']]
        depth = pd.concat([forward_depth, reverse_depth], axis=0, ignore_index=True)
        
        # Write output before sorting
        raw_output = out_path+'genomic_bin_count_binSize_'+str(windowSize)+'_stepSize_'+str(stepSize)+'_full_overlap.bed'
        depth.to_csv(raw_output, sep='\t', index=False, header=False)
        
        # Sort and cleanup
        sorted_output = out_path+'genomic_bin_count_binSize_'+str(windowSize)+'_stepSize_'+str(stepSize)+'_full_overlap.sorted.bed'
        safe_bedsort(raw_output, sorted_output)
        safe_remove(raw_output)
        
        # Cleanup temp unzipped alignment file
        if(os.path.exists(self.prefix + 'alignment.sorted.bed'+'.gz') and os.path.exists(self.prefix + 'alignment.sorted.bed')):
            safe_remove(self.prefix + 'alignment.sorted.bed')
            
        del windows, forward, reverse, forward_file, reverse_file, forward_depth, reverse_depth, depth
        