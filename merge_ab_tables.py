#!/usr/bin/env python3

import sys

"""
Merge script producing a final 14-column file with columns:
  1)  SNP_ID
  2)  CHR
  3)  POS
  4)  Sample_ID
  5)  VCF_REF
  6)  VCF_ALT
  7)  AlleleA
  8)  AlleleB
  9)  GenotypeAB
  10) Allele1_AB
  11) Allele2_AB
  12) Strand
  13) SurroundingSequence
  14) SourceGT

- file1_table.txt has 13 columns (no SurroundingSequence).
- file2_table.txt has 14 columns (including SurroundingSequence).

We insert "NA" for SurroundingSequence in file1 rows, ensuring final_merged.txt always has 14 columns.
"""

def merge_files(file1, file2, output):
    # Final header: 14 columns
    header = (
        "SNP_ID\tCHR\tPOS\tSample_ID\tVCF_REF\tVCF_ALT\t"
        "AlleleA\tAlleleB\tGenotypeAB\tAllele1_AB\tAllele2_AB\t"
        "Strand\tSurroundingSequence\tSourceGT\n"
    )

    with open(output, 'w') as fout:
        fout.write(header)

        # Process file1 (13 columns)
        with open(file1, 'r') as f1:
            f1.readline()  # skip file1 header
            for line in f1:
                line = line.strip()
                if not line:
                    continue
                cols = line.split('\t')
                if len(cols) < 13:
                    # skip malformed lines
                    continue
                # file1 columns => 13:
                # 0:SNP_ID,1:CHR,2:POS,3:Sample_ID,4:VCF_REF,5:VCF_ALT,
                # 6:AlleleA,7:AlleleB,8:GenotypeAB,9:Allele1_AB,10:Allele2_AB,
                # 11:Strand,12:SourceGT

                # We want to produce 14 columns by inserting "NA" for SurroundingSequence
                # at column index 12, just before SourceGT (index 13).
                # So let's do:
                new_cols = cols[:12] + ["NA"] + cols[12:]  # Now we have 14 columns
                fout.write("\t".join(new_cols) + "\n")

        # Process file2 (14 columns) if it has lines
        with open(file2, 'r') as f2:
            f2.readline()  # skip file2 header
            for line in f2:
                line = line.rstrip('\n')
                if not line.strip():
                    continue
                # file2 lines => 14 columns
                # We'll write them as-is
                fout.write(line + "\n")

def main():
    if len(sys.argv) < 4:
        print("Usage: python merge_ab_tables.py file1_table.txt file2_table.txt final_merged.txt", file=sys.stderr)
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output = sys.argv[3]
    merge_files(file1, file2, output)

if __name__ == "__main__":
    main()
