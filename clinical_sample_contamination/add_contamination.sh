# Mixes two clincial sample BAM files to simulate contamination at various levels
# Usage:  ./add_contamination.sh NA12878.bam NA19240.bam test


#!/bin/bash
set -euo pipefail

BASE="$1"      # base sample BAM (first argument)
CONTAM="$2"    # contaminant BAM (second argument)
PREFIX="$3"    # output prefix (third argument)

# List of contamination levels to simulate
ALPHAS=(0.01 0.02 0.03 0.04 0.05) 

for ALPHA in "${ALPHAS[@]}"; do
    echo " Simulating contamination alpha = $ALPHA"

    # Fractions for subsampling from the base BAM
    BASE_FRAC=$(python3 - <<EOF
alpha=float("$ALPHA")
frac = 1 - alpha
print(str(frac).split('.')[1])   # remove the leading "0."
EOF
)
    # Fractions for subsampling from the contaminant BAM
    CONTAM_FRAC=$(python3 - <<EOF
print("$ALPHA".split('.')[1])
EOF
)
    # Define output filenames
    CONTAM_LEVEL=$(printf "%03d" $(echo "$ALPHA*100" | bc))   # converts decimal to integer e.g. 0.03 → 003
    BASE_SUB="${PREFIX}.${CONTAM_LEVEL}.base.sub.bam"         # subsampled base BAM
    CONTAM_SUB="${PREFIX}.${CONTAM_LEVEL}.contam.sub.bam"     # subsampled contaminant BAM
    MERGED_RAW="${PREFIX}.${CONTAM_LEVEL}.raw.bam"            # merged unsorted BAM
    MERGED_SORTED="${PREFIX}.${CONTAM_LEVEL}.sorted.bam"      # merged sorted BAM
    FINAL="${PREFIX}.${CONTAM_LEVEL}.mixed.final.bam"         # final BAM with replaced read groups


    # Subsample base and contam BAMs at read level to give mixing fractions
    echo "→ Subsampling base ($BASE_FRAC)"
    samtools view -s 100."$BASE_FRAC" -b "$BASE" -o "$BASE_SUB"  # random seed 100

    echo "→ Subsampling contam ($CONTAM_FRAC)"
    samtools view -s 200."$CONTAM_FRAC" -b "$CONTAM" -o "$CONTAM_SUB"  # random seed 200

    echo "→ Merging"
    samtools merge -o "$MERGED_RAW" "$BASE_SUB" "$CONTAM_SUB"

    echo "→ Sorting"
    samtools sort "$MERGED_RAW" -o "$MERGED_SORTED"
    samtools index "$MERGED_SORTED"

    # Replace read groups in the merged BAM to have consistent RG info for compatibility with variant callers
    
    echo "Detecting platform (PL) from base BAM" # Detect sequencing platform from base BAM RG header
    PLATFORM=$(samtools view -H "$BASE" | grep "^@RG" | head -n 1 | sed 's/.*PL:\([^ \t]*\).*/\1/')

    echo "   Found platform: $PLATFORM"

    echo "Replacing read groups"
    java -jar picard.jar AddOrReplaceReadGroups \
        I="$MERGED_SORTED" \
        O="$FINAL" \
        RGID="SIM" \
        RGLB="SIMLIB" \
        RGPL="$PLATFORM" \
        RGPU="SIMUNIT" \
        RGSM="SIMSAMPLE"

    samtools index "$FINAL"

    echo "Finished alpha = $ALPHA → $FINAL"
    echo

    # Cleanup
    rm "$BASE_SUB" "$CONTAM_SUB" "$MERGED_RAW" "$MERGED_SORTED" "${MERGED_SORTED}.bai"
done
