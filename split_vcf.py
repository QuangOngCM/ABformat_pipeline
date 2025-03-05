#!/usr/bin/env python3

import sys

# Usage:
#   python split_vcf_step1.py input.vcf file1.vcf file2.vcf
#
# File 1: A↔C/G or T↔C/G
# File 2: A↔T or C↔G

def is_file1_snp(ref, alt):
    """
    Return True if SNP is A-C, A-G, T-C, or T-G (in either direction).
    """
    pair = {ref, alt}
    # Possible sets for File 1
    file1_sets = [
        {'A', 'C'},
        {'A', 'G'},
        {'T', 'C'},
        {'T', 'G'}
    ]
    return pair in file1_sets

def is_file2_snp(ref, alt):
    """
    Return True if SNP is A-T or C-G (in either direction).
    """
    pair = {ref, alt}
    # Possible sets for File 2
    file2_sets = [
        {'A', 'T'},
        {'C', 'G'}
    ]
    return pair in file2_sets

def main():
    if len(sys.argv) < 4:
        print("Usage: python split_vcf_step1.py input.vcf file1.vcf file2.vcf", file=sys.stderr)
        sys.exit(1)

    input_vcf = sys.argv[1]
    file1_vcf = sys.argv[2]
    file2_vcf = sys.argv[3]

    with open(input_vcf, 'r') as fin, \
         open(file1_vcf, 'w') as f1, \
         open(file2_vcf, 'w') as f2:

        for line in fin:
            # Pass header lines directly into both files
            if line.startswith("#"):
                f1.write(line)
                f2.write(line)
                continue

            # Parse the line
            cols = line.strip().split("\t")
            if len(cols) < 5:
                continue  # skip any malformed line

            ref = cols[3].upper()
            alt = cols[4].upper()

            # If there's more than one ALT allele (comma-separated),
            # skip or handle carefully. For simplicity, skip them here:
            if "," in alt:
                continue

            # Check if this is an A/C/G/T SNP
            # We'll also skip lines where REF or ALT is not a single base A/C/G/T.
            valid_bases = {'A','C','G','T'}
            if ref not in valid_bases or alt not in valid_bases:
                continue

            # Classify
            if is_file1_snp(ref, alt):
                f1.write(line)
            elif is_file2_snp(ref, alt):
                f2.write(line)
            else:
                # It's not in either category (e.g., could be other combos or multi-allelic).
                # We simply skip it here.
                pass

if __name__ == "__main__":
    main()
