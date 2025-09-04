"""
This script simulates contamination in sequencing data by mixing reads from a base BAM file
with reads from a contaminant BAM file at specified contamination levels (alphas). It outputs
new BAM files for each contamination level, along with a truth table CSV file documenting the
composition of each output BAM.
For more accurate contamination levels, BAM inputs should be filtered only for mapped reads, 
and have duplicates removed. BAM inputs must be indexed and sorted.

Usage: python simulate_contamination.py base_sample.bam contam_sample.bam
"""
#!/usr/bin/env python3
import pysam
import pysam.samtools
import random
import csv
import os

def get_aligned_bases(bamfile):
    """
    Return total aligned bases in a BAM using pysam.samtools.stats.
    """
    stats = pysam.samtools.stats(bamfile)
    for line in stats.splitlines():
        if line.startswith("SN") and "bases mapped (cigar)" in line:
            fields = line.strip().split()
            # Find the last numeric field in the line
            for field in reversed(fields):
                if field.isdigit():
                    return int(field)
    raise ValueError(f"Could not find aligned bases in {bamfile}")

def simulate_contamination(base_bam, contam_bam, 
                           alphas=[0.01, 0.02, 0.05, 0.1, 0.2],
                           seed=42, tag_contam=True, truth_csv="truth_table.csv"): # define function and inputs
    
    random.seed(seed) # random seed fixed for reproducibility in generating random samples

    # Generate output prefix from base and contam BAM names
    base_prefix = os.path.basename(base_bam).split(".")[0]
    contam_prefix = os.path.basename(contam_bam).split(".")[0]
    out_prefix = f"{base_prefix}_{contam_prefix}"

    # Open BAMs for reading (necessary for manipulating with pysam)
    base = pysam.AlignmentFile(base_bam, "rb") # rb = read binary (i.e. read bam)
    contam = pysam.AlignmentFile(contam_bam, "rb")

    # Count number of reads in each BAM for sampling calculations
    base_readcount = base.count(until_eof=True)
    contam_readcount = contam.count(until_eof=True)
    B1 = get_aligned_bases(base_bam) # total aligned bases in base BAM
    B2 = get_aligned_bases(contam_bam) # total aligned bases in contaminant BAM

    print(f"Base reads: {base_readcount}, Contaminant reads: {contam_readcount}")
    print(f"Base aligned bases: {B1}, Contaminant aligned bases: {B2}")

    # Prepare truth table file
    with open(truth_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["alpha", "n_base_reads", "n_contam_reads", "base_level_fraction", "read_level_fraction", "out_bam"])  # header row

        for alpha in alphas:
            base.reset() # move pointer back to start of BAM for re-iteration (avoids reopening file)
            contam.reset()

            # Probability of selecting individual reads from each BAM
            p_base = 1 - alpha
            p_contam = alpha * (B1 / B2)

            # Count of how many reads/bases from each BAM have been selected
            n_base_selected = 0
            n_contam_selected = 0
            total_base_bases = 0
            total_contam_bases = 0

            # Output BAM file
            out_bam = f"{out_prefix}_alpha{alpha:.2f}.bam"
            with pysam.AlignmentFile(out_bam, "wb", header=base.header) as out: # wb = write binary. imports header from base BAM
                # Sample base reads probabilistically
                for read in base.fetch(until_eof=True):
                    if random.random() < p_base:
                        out.write(read)
                        n_base_selected += 1
                        total_base_bases += read.reference_length

                # Sample contaminant reads probabilistically
                for read in contam.fetch(until_eof=True):
                    if random.random() < p_contam:
                        if tag_contam:
                            read.query_name = read.query_name + ":CONTAM"
                        out.write(read)
                        n_contam_selected += 1
                        total_contam_bases += read.reference_length

            # Calculate contamination fractions
            base_level_fraction = total_contam_bases / (total_base_bases + total_contam_bases) if (total_base_bases + total_contam_bases) > 0 else 0
            read_level_fraction = n_contam_selected / (n_base_selected + n_contam_selected) if (n_base_selected + n_contam_selected) > 0 else 0

            # Sort and index
            sorted_bam = out_bam.replace(".bam", ".sorted.bam")
            pysam.sort("-o", sorted_bam, out_bam, catch_stdout=False)
            pysam.index(sorted_bam)
            os.remove(out_bam)

            print(f"[α={alpha}] Selected {n_base_selected} base reads + {n_contam_selected} contaminant reads -> {sorted_bam}")

            writer.writerow([alpha, n_base_selected, n_contam_selected, round(base_level_fraction,5), round(read_level_fraction,5), sorted_bam]) # write to truth table

    base.close()
    contam.close()
    print(f"Truth table saved to {truth_csv}")


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage: python simulate_contamination.py base_sample.bam contam_sample.bam")
        sys.exit(1)

    base_bam = sys.argv[1]
    contam_bam = sys.argv[2]

    simulate_contamination(
        base_bam = base_bam,
        contam_bam = contam_bam)
    print("Contamination simulation complete.")

