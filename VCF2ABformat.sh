#!/usr/bin/env bash
#
# run_pipeline.sh
#
# Usage:
#   ./run_pipeline.sh <input_vcf> <reference_fasta>
#
# Example:
#   ./run_pipeline.sh test.vcf reference.fasta

set -e  # Exit immediately on error

# 1) Parse arguments
if [ $# -lt 2 ]; then
  echo "Usage: $0 <input_vcf> <reference_fasta>"
  exit 1
fi

INPUT_VCF="$1"
REFERENCE="$2"

# 2) Split the VCF into file1.vcf (non-ambiguous) and file2.vcf (ambiguous)
echo "==> Splitting $INPUT_VCF into file1.vcf and file2.vcf..."
python split_vcf.py "$INPUT_VCF" file1.vcf file2.vcf

# 3) Convert non-ambiguous SNPs (file1.vcf) => file1_table.txt
echo "==> Converting non-ambiguous SNPs (file1.vcf) to file1_table.txt..."
python nonambigous_convert.py file1.vcf file1_table.txt

# 4) Convert ambiguous SNPs (file2.vcf) => file2_table.txt (requires reference FASTA)
echo "==> Converting ambiguous SNPs (file2.vcf) to file2_table.txt..."
python ambigous_convert.py file2.vcf "$REFERENCE" file2_table.txt

# 5) Merge file1_table.txt + file2_table.txt => final_merged.txt
echo "==> Merging file1_table.txt and file2_table.txt into final_merged.txt..."
python merge_ab_tables.py file1_table.txt file2_table.txt final_merged.txt

# 6) Make genotype & map files from final_merged.txt
echo "==> Generating genotype_file.txt, snp_map.txt, snp_info_file.txt..."
python make_genotype_and_map.py final_merged.txt genotype_file.txt snp_map.txt snp_info_file.txt

# 7) (Optional) Make a FImpute-style genotype file
# Uncomment if you want to create genotype_file_fimpute.txt as well:
# echo "==> Creating genotype_file_fimpute.txt..."
# python make_genotype_file_fimpute.py genotype_file.txt snp_info_file.txt genotype_file_fimpute.txt

echo "==> All done!"
