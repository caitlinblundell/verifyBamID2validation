# Script to filter, sort, and index sample files for in silico contamination, then downsample to 12M aligned reads
# Using samtools CLI to preserve read group information
# Usage: ./prepare_bams.sh file1.bam file2.bam ...

#!/bin/bash
set -euo pipefail

# Set target number of aligned reads after downsampling
TARGET_ALIGNED_READS=12000000  

for bam in "$@"; do  # Iterate over all input BAM files
    # Extract sample info from BAM filename
    fname=$(basename "$bam")  # Get filename without path
    sample="${fname%.bam}"    # Remove .bam extension

    filtered="${sample}.filtered.no_dups.bam"
    downsampled="${sample}.filtered.no_dups.12M.bam"
    sorted="${sample}.filtered.no_dups.12M.sorted.bam"

    echo "Processing sample: $sample"

    # Filter out unmapped/duplicate/secondary/supplementary reads
    echo "Filtering $bam"
    samtools view -b -h -F 3332 "$bam" -o "$filtered"
    
    # Count mapped reads for downsampling fraction
    MAPPED_READS=$(samtools view -c -F 4 "$filtered")
    echo "   Mapped reads before downsampling: $MAPPED_READS"

    # Calculate downsampling fraction
    if (( MAPPED_READS <= TARGET_ALIGNED_READS )); then
        echo "Sample has fewer than target reads; skipping downsampling."
        DOWNSAMPLE_FRAC=1.0
    else
        DOWNSAMPLE_FRAC=$(python3 - <<EOF
frac = $TARGET_ALIGNED_READS / $MAPPED_READS
# remove leading "0." for samtools view -s option
print(str(frac).split(".")[1])
EOF
)
        echo "   Downsampling fraction: 0.$DOWNSAMPLE_FRAC"
    fi

    # Downsample to target aligned reads:
    echo "Downsampling to $TARGET_ALIGNED_READS aligned reads"
    samtools view -s "42.$DOWNSAMPLE_FRAC" -b "$filtered" -o "$downsampled" #random seed 42

    # Sort + index
    echo "Sorting"
    samtools sort "$downsampled" -o "$sorted"

    echo "Indexing"
    samtools index "$sorted"

    # Clean up
    echo "Removing intermediate file"
    rm "$filtered" "$downsampled"

done
