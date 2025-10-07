#!/bin/bash
set -e

# This script merges all 1KGP individual chr VCF files into a single VCF file using bcftools.

VCF_DIR="1kgp_phase3_vcfs"
PANEL_BED="Pan5272_data.bed"
SUBSET_VCF_DIR="VCF_subset"

mkdir -p $SUBSET_VCF_DIR

# Log file
VARIANT_LOG_FILE="variant_counts_per_chr.txt"
echo "Chromosome   NumVariants" > $VARIANT_LOG_FILE

# Subset each with CP2 .bed and then compress/index
for vcf in $VCF_DIR/ALL.chr*.vcf.gz; do
    vcf_prefix=$(basename $vcf | cut -d. -f2)  # chr1, chr2, etc.
    echo "Processing $vcf_prefix"
    bcftools view -R $PANEL_BED $vcf -Oz -o $SUBSET_VCF_DIR/$vcf_prefix.subset.vcf.gz
    bcftools index $SUBSET_VCF_DIR/$vcf_prefix.subset.vcf.gz

    #logging
    num_variants=$(bcftools view -H $SUBSET_VCF_DIR/$vcf_prefix.subset.vcf.gz | wc -l)
    echo -e "$vcf_prefix\t$num_variants" >> $VARIANT_LOG_FILE
done

# Merge all subset VCFs into a single VCF for filtering with bcftools
bcftools concat -Oz -o merged_1kgp_panel.vcf $SUBSET_VCF_DIR/*.subset.vcf.gz
bcftools index merged_1kgp_panel.vcf.gz