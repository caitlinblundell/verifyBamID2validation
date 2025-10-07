#!/bin/bash
set -e

VCF_PANEL="merged_1kgp_panel.vcf.gz"

# Step 1: Retain only biallelic SNPs with genotype information and filter=PASS
bcftools view -m2 -M2 -v snps -f PASS -i 'FMT/GT!="."' $VCF_PANEL -Oz -o step1_filtered_snps.vcf.gz

# Step 2: Keep only variants with genotype missing rate below 20%
bcftools +fill-tags step1_filtered_snps.vcf.gz -Oz -o step2_filled_tags.vcf.gz -- -t F_MISSING
bcftools view -i 'F_MISSING<0.2' step2_filled_tags.vcf.gz -Oz -o step2_filtered_1kgp_panel.vcf.gz

# Count  variants
current_variants=$(bcftools view -H step2_filtered_1kgp_panel.vcf.gz | wc -l)
echo "✅ Total variants pre-masking and filtering for MAF: $current_variants"

MAF_CUTOFF=0.005

# Step 3: Mask genotypes according to 1KGP strict mask
bedtools intersect -header -v -a step2_filtered_1kgp_panel.vcf.gz -b 20141020.strict_mask.whole_genome.bed > step3_remove_masked.vcf
bgzip step3_unmasked.vcf
bcftools index step3_unmasked.vcf.gz

# Step 4: Filter variants based on MAF
bcftools +fill-tags step3_remove_masked.vcf.gz -Oz -o step4_filter_MAF.vcf.gz -- -t MAF
bcftools view -i 'MAF>0.005' step4_filter_MAF.vcf.gz -Oz -o 1KGP_CP2_panel.vcf.gz

# Index final VCF
bcftools index 1KGP_CP2_panel.vcf.gz

# Count total variants
total_variants=$(bcftools view -H 1KGP_CP2_panel.vcf.gz | wc -l)
echo "✅ Total variants after all filters: $total_variants"