#!/bin/bash

# script was run interactively - but could ideally be included in the workflow

source activate rnaseq
cd data/external/strains/

#### VCF files ####

# Download VCFs (and their indexes) - SNPs and indels for CAST/EiJ strain
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz.tbi
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/CAST_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/CAST_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz.tbi

# concatenate SNPs and indels (the input files must be sorted by chr and position)
# retain only variants that pass quality filters
# and turn individual to a "heterozygote" (F1 hybrid) with sed substitutions
bcftools concat \
  -a --remove-duplicates -Ov \
  CAST_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz |\
  bcftools filter --include 'FILTER="PASS"' -Ov |\
  sed 's/\t1\/1:/\t0\/1:/' |\
  sed 's/CAST_EiJ$/B6xCAST/' |\
  bgzip > B6xCAST.snps_and_indels.vcf.gz

bcftools index B6xCAST.snps_and_indels.vcf.gz

# create BED file (for downstream analysis in R, etc.)
# sorting by chromosome and position (numeric)
# merging overlapping SNPs and Indels (the -d -1 is so that adjoining SNPs don't get merged)
bcftools query \
  --format '%CHROM\t%POS0\t%END\t%ID-%TYPE\n' \
  B6xCAST.snps_and_indels.vcf.gz | \
  sort -k1,1 -k2n,2n | \
  bedtools merge -c 4 -o collapse -d -1 | \
  gzip > B6xCAST.snps_and_indels.merged.bed.gz

