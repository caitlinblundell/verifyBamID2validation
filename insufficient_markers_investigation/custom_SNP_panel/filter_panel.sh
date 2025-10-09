#!/bin/bash
set -e

# This script filters the merged 1kgp panel according to the criteria specified by VBID2, producing vbid2_1KGP_CP2_panel.vcf.gz. 
# It also filters based on a MAF cutoff in case this more stringent panel should be used instead

VCF_PANEL="1KGP_CP2_intersect/merged_1kgp_panel.vcf.gz"
TEMP_FILTER_STEPS_DIR="interim_filter_files"

mkdir -p $TEMP_FILTER_STEPS_DIR

# Step 1: Retain only biallelic SNPs with genotype information and filter=PASS
bcftools view -m2 -M2 -v snps -f PASS -i 'FMT/GT!="."' $VCF_PANEL -Oz -o $TEMP_FILTER_STEPS_DIR/step1_filtered_snps.vcf.gz

# Step 2: Keep only variants with genotype missing rate below 20%
bcftools +fill-tags $TEMP_FILTER_STEPS_DIR/step1_filtered_snps.vcf.gz -Oz -o $TEMP_FILTER_STEPS_DIR/GT_filled_tags.vcf.gz -- -t F_MISSING
bcftools view -i 'F_MISSING<0.2' $TEMP_FILTER_STEPS_DIR/GT_filled_tags.vcf.gz -Oz -o $TEMP_FILTER_STEPS_DIR/step2_filtered_1kgp_panel.vcf.gz

# Count  variants
step2_variants=$(bcftools view -H $TEMP_FILTER_STEPS_DIR/step2_filtered_1kgp_panel.vcf.gz | wc -l)
echo "Total variants pre-masking and filtering for MAF: $step2_variants"

# Step 3: Mask genotypes according to 1KGP strict mask
bedtools intersect -header -v -a $TEMP_FILTER_STEPS_DIR/step2_filtered_1kgp_panel.vcf.gz -b 1KGP_mask/20141020.strict_mask.whole_genome.bed > vbid2_1KGP_CP2_panel.vcf
bgzip vbid2_1KGP_CP2_panel.vcf
bcftools index vbid2_1KGP_CP2_panel.vcf.gz

# Count  variants
step3_variants=$(bcftools view -H vbid2_1KGP_CP2_panel.vcf.gz | wc -l)
echo "Total variants after masking: $step3_variants"

# Step 4: Filter variants based on MAF
MAF_CUTOFF=0.005
bcftools +fill-tags vbid2_1KGP_CP2_panel.vcf.gz -Oz -o $TEMP_FILTER_STEPS_DIR/step4_filter_MAF.vcf.gz -- -t MAF
bcftools view -i "MAF>${MAF_CUTOFF}" $TEMP_FILTER_STEPS_DIR/step4_filter_MAF.vcf.gz -Oz -o MAF_filter_1KGP_CP2_panel.vcf.gz
bcftools index MAF_filter_1KGP_CP2_panel.vcf.gz

# Count total variants
total_variants=$(bcftools view -H MAF_filter_1KGP_CP2_panel.vcf.gz | wc -l)
echo "Total variants after MAF filter: $total_variants"