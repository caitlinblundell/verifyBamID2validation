#!/bin/bash

# Index file
INDEX="20130502.phase3.low_coverage.alignment.index"

# Loop over all BAM files in the current directory
for bamfile in *.bam; do
    # Check if the file exists
    if [[ ! -f "$bamfile" ]]; then
        continue
    fi

    # Extract the corresponding MD5 from the index
    # The index has server paths, so match by basename
    bam_md5=$(awk -v fname="$bamfile" 'NR>1 && $1 ~ fname {print $2}' "$INDEX")

    if [[ -z "$bam_md5" ]]; then
        echo "  WARNING: $bamfile not found in index"
        continue
    fi

    # Check MD5
    echo "Checking $bamfile ..."
    md5sum "$bamfile" | awk -v expected="$bam_md5" '{
        if ($1 == expected) {
            print "  OK:", $2
        } else {
            print "  MISMATCH:", $2, "(expected " expected ", got " $1 ")"
        }
    }'
done
