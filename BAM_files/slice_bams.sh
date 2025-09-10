#!/bin/bash

# Slice all BAMs in the current directory according to the vbid2 BED file
# Usage: ./slice_bams.sh

BEDFILE="1000g.phase3.100k.b37.vcf.gz.dat.bed"

# Loop over all BAM files in the current directory
for bam in *.bam; do
    # Skip if no BAM files found
    [[ -f "$bam" ]] || continue

    # Prepend "slice." to the original filename
    out="slice.$bam"

    echo "Slicing $bam using $BEDFILE → $out ..."

    # Slice BAM with samtools, keep header
    samtools view -b -h -L "$BEDFILE" "$bam" -o "$out"

    # index the sliced BAM
    samtools index "$out"

    echo "Done: $out"
done