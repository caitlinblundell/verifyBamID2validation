#!/bin/bash
set -e

VCF_PANEL="merged_1kgp_panel.vcf.gz"

# Step 1: Retain only biallelic SNPs
bcftools view -m2 -M2 -v snps $VCF_PANEL -Oz -o step1_biallelic_snps.vcf.gz

# Step2 2: Keep only filter=PASS variants
bcftools view -f PASS step1_biallelic_snps.vcf.gz -Oz -o step2_filter_pass.vcf

#Step 3: Keep only variants with GT field
bcftools view -i 'FMT/GT!="."' step2_filter_pass.vcf -Oz -o step3_with_GT.vcf.gz

# Step 4: Keep only variants with genotype missing rate below 20%
bcftools +fill-tags step3_with_GT.vcf.gz -- -t F_MISSING -Oz -o step4_filled_tags.vcf.gz
bcftools view -i 'F_MISSING<0.2' step4_filled_tags.vcf.gz -Oz -o filtered_1kgp_panel.vcf.gz

# Index final filtered VCF
bcftools index filtered_1kgp_panel.vcf.gz