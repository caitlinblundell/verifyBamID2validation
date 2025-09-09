"""
This script prepares downloaded BAM files for in silico contamination by filtering out unmapped and duplicate reads. 
The output BAMs are sorted and indexed, and named accordingly.

The input BAM filenames must follow the 1000 Genomes Project low coverage WGS naming convention:
<sampleID>.mapped.ILLUMINA.bwa.<ancestry>.{other_info}.bam

Usage: python prepare_bams.py sample1.bam sample2.bam
"""

#!/usr/bin/env python3
import pysam 
import pysam.samtools
import os

def extract_sample_info(bamfile):
    """
    Extract sample name and ancestry from BAM file header.
    """
    filename = os.path.basename(bamfile)
    parts = filename.split(".") # split filename by "." as per 1KGP naming convention
    sample_id = parts[0]  # first part is sample ID
    ancestry = parts[4]  # fifth part is ancestry code
    if len(ancestry) != 3:
        raise ValueError(f"Unexpected filename format for {bamfile}, cannot extract ancestry.")
    return sample_id, ancestry

def filter_bam(input_bams):
    """Filter unmapped/duplicate reads, sort, index, and rename BAM."""
    for input_bam in input_bams:
            sample_id, ancestry = extract_sample_info(input_bam)  # get ID and ancestry
            filtered_bam = f"{sample_id}.filtered_mapped.no_dups.{ancestry}.bam"
            sorted_bam = f"{sample_id}.filtered_mapped.no_dups.sorted.{ancestry}.bam"

            # filter out unmapped and duplicate reads
            pysam.samtools.view("-b", "-h", "-F", "1028", input_bam, "-o", filtered_bam, catch_stdout=False)
            print(f"{input_bam} filtered -> {filtered_bam}")
            # sort and index
            pysam.sort(filtered_bam, "-o", sorted_bam, catch_stdout=False)
            pysam.index(sorted_bam)
            print(f"{filtered_bam} sorted and indexed.")

            # Remove intermediate filtered BAM
            os.remove(filtered_bam)

            # Remove original BAM and its index if they exist
            bai_file = input_bam + ".bai"
            if os.path.exists(sorted_bam):  # only delete originals if sorted BAM exists
                if os.path.exists(input_bam):
                    os.remove(input_bam)
                if os.path.exists(bai_file):
                    os.remove(bai_file)
                print(f"Original BAM {input_bam} and index {bai_file} removed.\n")