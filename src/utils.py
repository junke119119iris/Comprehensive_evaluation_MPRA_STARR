# === Standard Library ===
import os
from subprocess import Popen, PIPE, STDOUT, run

# === Bioinformatics Tools ===
import pybedtools

def gzip_file(file):
    """
    Compress a file using gzip.
    
    Args:
        file (str): Path to the file to be compressed.
    """
    
    os.system("gzip {}".format(file))

def gunzip_file(file_gz, file=None, keep=True):
    
    """
    Decompress a gzip file.

    Args:
        file_gz (str): Path to the gzipped file.
        file (str, optional): Output path for decompressed content. If not provided, decompresses in place.
        keep (bool): If True and output file is specified, preserves the original .gz file.
    """
    if keep and file:
        os.system(f"gunzip < {file_gz} > {file}")
    else:
        os.system(f"gunzip {file_gz}")

def safe_bedsort(input_file, output_file):
    """
    Sort a BED file by chromosome and start coordinate.

    Args:
        input_file (str): Path to input BED file.
        output_file (str): Path to write the sorted BED file.
    """
    
    os.system(f"sort -k1,1 -k2,2n {input_file} > {output_file}")

def safe_remove(file):
    """
    Remove a file if it exists.

    Args:
        file (str): Path to the file to be removed.
    """
    
    if os.path.exists(file):
        os.remove(file)

def set_pybedtools_tmp_dir(tmp_dir):
    """
    Set a temporary directory for pybedtools.

    Args:
        tmp_dir (str): Path to temporary directory.
    """
    
    pybedtools.helpers.set_tempdir(tmp_dir)

def set_dir(path):
    """
    Create a directory if it doesn't already exist.

    Args:
        path (str): Directory path to create.
    """
    
    if not os.path.exists(path):
        os.mkdir(path)

def get_row_number(file):
    """
    Count the number of lines (rows) in a file.

    Args:
        file (str): Path to the file.

    Returns:
        int: Number of lines in the file.
    """
    
    results = run(['wc', '-l', file], stdout=PIPE, stderr=PIPE, universal_newlines=True)
    return int(results.stdout.split()[0])

def get_row_number_zipped(file):
    
    """
    Count the number of lines in a file, including gzipped files.

    Args:
        file (str): Path to the (possibly gzipped) file.

    Returns:
        int: Number of lines in the file.
    """
    
    cmd = f"zcat {file} | wc -l" if file.endswith('.gz') else f"wc -l {file}"
    ps = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    output = ps.communicate()[0].decode("utf-8")
    return int(output.strip().split()[0])
